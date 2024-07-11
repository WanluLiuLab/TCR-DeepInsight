from collections import Counter
import torch
from torch import nn
from torch.nn import functional as F
import json
import tqdm
from abc import ABC
import datasets
from typing import Optional, Union, Tuple, Mapping
from transformers import PreTrainedModel
import numpy as np
from pathlib import Path
from tokenizers import Tokenizer
import pandas as pd

from scatlasvae.utils._tensor_utils import one_hot

from ._collator import TRabCollatorForVJCDR3, AminoAcidsCollator
from ._model import TRabModelingBertForPseudoSequence, TRabModelingBertForVJCDR3
from .._model_utils import EarlyStopping

class TrainerBase(ABC):
    def attach_train_dataset(self, train_dataset: Optional[datasets.Dataset]):
        self.train_dataset = train_dataset

    def attach_test_dataset(self, test_dataset: Optional[datasets.Dataset]):
        self.test_dataset = test_dataset

class TRabModelingBertForVJCDR3Trainer(TrainerBase):
    def __init__(
        self,
        model: TRabModelingBertForVJCDR3,
        collator: TRabCollatorForVJCDR3,
        train_dataset: Optional[datasets.Dataset] = None,
        test_dataset: Optional[datasets.Dataset] = None,
        optimizers: Tuple[torch.optim.Optimizer, torch.optim.lr_scheduler.LambdaLR] = (None, None),
        loss_weight: Mapping[str, float] = {"reconstruction_loss": 1},
        device: str = 'cuda'
    ) -> None:
        self.model = model
        self.optimizer, self.scheduler = optimizers
        self.device = device
        if self.device != model.device:
            self.device = model.device
        self.collator = collator 
        self.loss_weight = loss_weight
        self.attach_train_dataset(train_dataset)
        self.attach_test_dataset(test_dataset)

    def fit(
        self,
        *,
        max_epoch: int,
        max_train_sequence: int = 0,
        n_per_batch: int = 10,
        shuffle: bool = False,
        balance_label: bool = False,
        label_weight: Optional[torch.Tensor] = None,
        show_progress: bool = False,
        early_stopping: bool = False
    ):
        """Main training entry point"""
        loss = []
        self.model.train()
        if early_stopping:
            early_stopper = EarlyStopping()
        if max_train_sequence == 0:
            max_train_sequence = len(self.train_dataset)
        if shuffle: self.train_dataset.shuffle()
        for i in range(1, max_epoch + 1):
            epoch_total_loss = []
            if show_progress:
                pbar = tqdm.tqdm(total=max_train_sequence // n_per_batch)

            for j in range(0, max_train_sequence, n_per_batch):
                self.optimizer.zero_grad()
                epoch_indices = np.array(
                    self.train_dataset[j:j+n_per_batch]['input_ids']
                )

                epoch_attention_mask = np.array(
                    self.train_dataset[j:j+n_per_batch]['attention_mask']
                )

                if 'cdr3_mr_mask' in self.train_dataset.features.keys():
                    cdr3_mr_mask = np.array(
                        self.train_dataset[j:j+n_per_batch]['cdr3_mr_mask']
                    )
                else:
                    cdr3_mr_mask = self.train_dataset[j:j+n_per_batch]['attention_mask']
                    
                epoch_token_type_ids = torch.tensor(
                    self.train_dataset[j:j+n_per_batch]['token_type_ids'], dtype=torch.int64
                ).to(self.device)

                epoch_indices_mlm, epoch_attention_mask_mlm = self.collator(
                    epoch_indices, epoch_attention_mask, cdr3_mr_mask
                )
                # print(tokenizer.convert_ids_to_tokens(torch.tensor(epoch_indices_mlm)))
                epoch_indices_mlm = epoch_indices_mlm.to(self.device)
                epoch_attention_mask_mlm = epoch_attention_mask_mlm.to(self.device)
                epoch_label_ids = None

                if 'tcr_label' in self.train_dataset.features.keys():
                    epoch_label_ids = torch.tensor(
                        self.train_dataset[j:j+n_per_batch]['tcr_label'], dtype=torch.int64
                    ).to(self.device)
                    if balance_label:
                        label_counter = Counter(epoch_label_ids.cpu().numpy())
                        min_label_number = np.min(list(label_counter.values()))
                        min_label = {v:k for k,v in label_counter.items()}[min_label_number]
                        min_indices = np.argwhere(epoch_label_ids.cpu().numpy() == min_label).squeeze()
                        for k in label_counter.keys():
                            if k != min_label:
                                indices = np.argwhere(epoch_label_ids.cpu().numpy() == k).squeeze()
                                np.random.shuffle(indices)
                                indices = indices[:min_label_number]
                                min_indices = np.concatenate([min_indices, indices])
                        epoch_indices_mlm = epoch_indices_mlm[min_indices]
                        epoch_attention_mask_mlm = epoch_attention_mask_mlm[min_indices]
                        epoch_label_ids = epoch_label_ids[min_indices]
                        epoch_token_type_ids = epoch_token_type_ids[min_indices]
                        epoch_indices = epoch_indices[min_indices]


                self.optimizer.zero_grad()
                output = self.model.forward(
                    input_ids = epoch_indices_mlm,
                    attention_mask = epoch_attention_mask_mlm,
                    token_type_ids =  epoch_token_type_ids,
                    labels = torch.tensor(epoch_indices, 
                    dtype=torch.int64).to(self.device)
                )

                if epoch_label_ids is not None:
                    prediction_loss = nn.BCEWithLogitsLoss(
                        weight=label_weight.to(self.device)
                    )(
                        output["prediction_out"], 
                        one_hot(epoch_label_ids.unsqueeze(1), self.model.labels_number)
                    )
                else: 
                    prediction_loss = torch.tensor(0.)

                total_loss = output["output"].loss * self.loss_weight['reconstruction_loss'] + prediction_loss
                
                self.optimizer.zero_grad()
                total_loss.backward()
                self.optimizer.step()
                epoch_total_loss.append( total_loss.item() )

                if show_progress:
                    pbar.update(1)
                    pbar.set_postfix({
                        f'Learning rate': self.optimizer.param_groups[0]['lr'],
                        f'Loss': np.mean(epoch_total_loss[-4:]),
                        f'Prediction': prediction_loss.item()
                    })

            if self.scheduler is not None:
                self.scheduler.step(np.mean(epoch_total_loss))

            if early_stopping:
                early_stopper(np.mean(epoch_total_loss))
                if early_stopper.early_stop:
                    print("Early stopping at epoch {}".format(i))
                    break
            if show_progress:
                pbar.close()

            print("epoch {} | total loss {:.2e}".format(i, np.mean(epoch_total_loss)))
            
            loss.append(np.mean(epoch_total_loss))
        self.model.train(False)
        return loss

    def evaluate(
        self,
        *,
        n_per_batch: int = 10,
        max_train_sequence: int = None,
        max_test_sequence: int = None,
        show_progress: bool = False
    ):
        """Main training entry point"""
        self.model.eval()
        ### Train
        
        all_train_result = {
            "aa": [],
            "aa_pred": [],
            "aa_gt": [],
            "av": [],
            "bv": []
        }
        if self.train_dataset is not None:
            max_train_sequence = len(self.train_dataset) if max_train_sequence is None else max_train_sequence
            if show_progress:
                pbar = tqdm.tqdm(total=max_train_sequence // n_per_batch)
            for j in range(0, max_train_sequence, n_per_batch):
                self.optimizer.zero_grad()
                epoch_indices = np.array(
                        self.train_dataset[j:j+n_per_batch]['input_ids']
                )

                epoch_attention_mask = np.array(
                        self.train_dataset[j:j+n_per_batch]['attention_mask']
                )

                if 'cdr3_mr_mask' in self.train_dataset.features.keys():
                    cdr3_mr_mask = np.array(
                        self.train_dataset[j:j+n_per_batch]['cdr3_mr_mask']
                    )
                else:
                    cdr3_mr_mask = self.train_dataset[j:j+n_per_batch]['attention_mask']
                        
                epoch_token_type_ids = torch.tensor(
                        self.train_dataset[j:j+n_per_batch]['token_type_ids'], dtype=torch.int64
                ).to(self.device)

                epoch_indices_mlm, epoch_attention_mask_mlm = self.collator(
                        epoch_indices, epoch_attention_mask, cdr3_mr_mask
                )
                # print(tokenizer.convert_ids_to_tokens(torch.tensor(epoch_indices_mlm)))
                epoch_indices_mlm = epoch_indices_mlm.to(self.device)
                epoch_attention_mask_mlm = epoch_attention_mask_mlm.to(self.device)
                epoch_label_ids = None
                if 'tcr_label' in self.train_dataset.features.keys():
                    epoch_label_ids = torch.tensor(
                                self.train_dataset[j:j+n_per_batch]['tcr_label'], dtype=torch.int64
                    ).to(self.device)


                self.optimizer.zero_grad()
                output = self.model.forward(
                            input_ids = epoch_indices_mlm,
                            attention_mask = epoch_attention_mask_mlm,
                            token_type_ids =  epoch_token_type_ids,
                            # tcr_label_ids = epoch_label_ids,
                            labels = torch.tensor(epoch_indices, 
                            dtype=torch.int64).to(self.device)
                )

                predictions = output['output']['logits'].topk(1)[1].squeeze()
                evaluate_mask = epoch_indices_mlm == 22
                ground_truth = torch.tensor(epoch_indices, dtype=torch.int64).to(self.device)
                evaluate_mask_trav = torch.zeros_like(evaluate_mask, dtype=torch.bool)
                evaluate_mask_trav[:, 0] = True
                evaluate_mask_trbv = torch.zeros_like(evaluate_mask, dtype=torch.bool)
                evaluate_mask_trbv[:, 48] = True
                evaluate_mask[:, [0, 48]] = False
                result = (predictions == ground_truth)[evaluate_mask].detach().cpu().numpy()
                all_train_result['aa'].append(result)
                all_train_result['aa_pred'].append(predictions[evaluate_mask].detach().cpu().numpy())
                all_train_result['aa_gt'].append(ground_truth[evaluate_mask].detach().cpu().numpy())
                result = (predictions == ground_truth)[evaluate_mask_trav].detach().cpu().numpy()
                all_train_result['av'].append(result)
                result = (predictions == ground_truth)[evaluate_mask_trbv].detach().cpu().numpy()
                all_train_result['bv'].append(result)
                if show_progress:
                    pbar.update(1)
            if show_progress:
                pbar.close()

        ### Test
        
        all_test_result = {
            "aa": [],
            "av": [],
            "bv": []
        }
        if self.test_dataset is not None:
            if show_progress:
                pbar = tqdm.tqdm(total=max_test_sequence // n_per_batch)
            max_test_sequence = len(self.test_dataset) if max_test_sequence is None else max_test_sequence
            for j in range(0, max_test_sequence, n_per_batch):
                self.optimizer.zero_grad()
                epoch_indices = np.array(
                        self.test_dataset[j:j+n_per_batch]['input_ids']
                )

                epoch_attention_mask = np.array(
                        self.test_dataset[j:j+n_per_batch]['attention_mask']
                )

                if 'cdr3_mr_mask' in self.test_dataset.features.keys():
                    cdr3_mr_mask = np.array(
                        self.test_dataset[j:j+n_per_batch]['cdr3_mr_mask']
                    )
                else:
                    cdr3_mr_mask = self.test_dataset[j:j+n_per_batch]['attention_mask']
                        
                epoch_token_type_ids = torch.tensor(
                        self.test_dataset[j:j+n_per_batch]['token_type_ids'], dtype=torch.int64
                ).to(self.device)

                epoch_indices_mlm, epoch_attention_mask_mlm = self.collator(
                        epoch_indices, epoch_attention_mask, cdr3_mr_mask
                )
                # print(tokenizer.convert_ids_to_tokens(torch.tensor(epoch_indices_mlm)))
                epoch_indices_mlm = epoch_indices_mlm.to(self.device)
                epoch_attention_mask_mlm = epoch_attention_mask_mlm.to(self.device)
                epoch_label_ids = None
                if 'tcr_label' in self.test_dataset.features.keys():
                    epoch_label_ids = torch.tensor(
                                self.test_dataset[j:j+n_per_batch]['tcr_label'], dtype=torch.int64
                    ).to(self.device)


                self.optimizer.zero_grad()
                output = self.model.forward(
                            input_ids = epoch_indices_mlm,
                            attention_mask = epoch_attention_mask_mlm,
                            token_type_ids =  epoch_token_type_ids,
                            # tcr_label_ids = epoch_label_ids,
                            labels = torch.tensor(epoch_indices, 
                            dtype=torch.int64).to(self.device)
                )

                predictions = output['output']['logits'].topk(1)[1].squeeze()
                evaluate_mask = epoch_indices_mlm == 22
                ground_truth = torch.tensor(epoch_indices, dtype=torch.int64).to(self.device)
                evaluate_mask_trav = torch.zeros_like(evaluate_mask, dtype=torch.bool)
                evaluate_mask_trav[:, 0] = True
                evaluate_mask_trbv = torch.zeros_like(evaluate_mask, dtype=torch.bool)
                evaluate_mask_trbv[:, 48] = True
                evaluate_mask[:, [0, 48]] = False
                result = (predictions == ground_truth)[evaluate_mask].detach().cpu().numpy()
                all_test_result['aa'].append(result)
                result = (predictions == ground_truth)[evaluate_mask_trav].detach().cpu().numpy()
                all_test_result['av'].append(result)
                result = (predictions == ground_truth)[evaluate_mask_trbv].detach().cpu().numpy()
                all_test_result['bv'].append(result)
                if show_progress:
                    pbar.update(1)
        
            if show_progress:
                pbar.close()
        return all_train_result, all_test_result
    
class TRabModelingBertForPseudoSequenceTrainer(TrainerBase):
    def __init__(
        self,
        model: TRabModelingBertForPseudoSequence,
        collator: AminoAcidsCollator,
        train_dataset: Optional[Union[datasets.Dataset, pd.DataFrame]] = None,
        test_dataset: Optional[Union[datasets.Dataset, pd.DataFrame]] = None,
        optimizers: Tuple[torch.optim.Optimizer, torch.optim.lr_scheduler.LambdaLR] = (None, None),
        loss_weight: Mapping[str, float] = {"reconstruction_loss": 1},
        device: str = 'cuda'
    ) -> None:
        self.model = model
        self.optimizer, self.scheduler = optimizers
        self.device = device
        if self.device != model.device:
            self.device = model.device
        self.collator = collator 
        self.loss_weight = loss_weight

        self.attach_train_dataset(train_dataset)
        self.attach_test_dataset(test_dataset)

    def fit(
        self,
        *,
        max_epoch: int,
        n_per_batch: int = 10,
        shuffle: bool = False,
        show_progress: bool = False,
        early_stopping: bool = False,
        cache_file: Optional[Path] = None
    ):
        loss = []
        self.model.train()
        if early_stopping:
            early_stopper = EarlyStopping()

        max_train_sequence = len(self.train_dataset)
        if shuffle: 
            self.train_dataset.shuffle()

        for i in range(1, max_epoch + 1):
            epoch_total_loss = []
            if show_progress:
                pbar = tqdm.tqdm(total=max_train_sequence // n_per_batch)

            for j in range(0, max_train_sequence, n_per_batch):

                self.optimizer.zero_grad()

                epoch_indices = np.array(
                    self.train_dataset[j:j+n_per_batch]['input_ids']
                )

                epoch_attention_mask = np.array(
                    self.train_dataset[j:j+n_per_batch]['attention_mask']
                )

                epoch_indices_mlm, epoch_attention_mask_mlm = self.collator(
                    epoch_indices, epoch_attention_mask
                )
                

                # print(tokenizer.convert_ids_to_tokens(torch.tensor(epoch_indices_mlm)))
                epoch_indices_mlm = epoch_indices_mlm.to(self.device)
                epoch_attention_mask_mlm = epoch_attention_mask_mlm.to(self.device)
                epoch_label_ids = None

                if 'tcr_label' in self.train_dataset.features.keys():
                    epoch_label_ids = torch.tensor(
                        self.train_dataset[j:j+n_per_batch]['tcr_label'], dtype=torch.int64
                    ).to(self.device)

                self.optimizer.zero_grad()
                output = self.model.forward(
                    input_ids = epoch_indices_mlm,
                    attention_mask = epoch_attention_mask_mlm,
                    labels = torch.tensor(
                        epoch_indices, 
                        dtype=torch.int64
                    ).to(self.device)
                )

                if epoch_label_ids is not None:
                    prediction_loss = nn.BCEWithLogitsLoss()(
                        output["prediction_out"], 
                        one_hot(epoch_label_ids.unsqueeze(1), self.model.labels_number)
                    )
                else: 
                    prediction_loss = torch.tensor(0.)

                total_loss = output["output"].loss * self.loss_weight['reconstruction_loss'] + prediction_loss
                
                self.optimizer.zero_grad()
                total_loss.backward()
                self.optimizer.step()
                epoch_total_loss.append( total_loss.item() )

                if show_progress:
                    pbar.update(1)
                    pbar.set_postfix({
                        f'Learning rate': self.optimizer.param_groups[0]['lr'],
                        f'Loss': np.mean(epoch_total_loss[-4:]),
                        f'Prediction': prediction_loss.item()
                    })

            if self.scheduler is not None:
                self.scheduler.step(np.mean(epoch_total_loss))

            if early_stopping:
                early_stopper(np.mean(epoch_total_loss))
                if early_stopper.early_stop:
                    print("Early stopping at epoch {}".format(i))
                    break
            if show_progress:
                pbar.close()

            print("epoch {} | total loss {:.2e}".format(i, np.mean(epoch_total_loss)))
            
            loss.append(np.mean(epoch_total_loss))
        self.model.train(False)
        return loss
    
    def evaluate(
        self, *,
        n_per_batch: int = 10,
        max_train_sequence: int = None,
        max_test_sequence: int = None,
        show_progress: bool = False
    ):
        self.model.eval()

        ### Train
        
        all_train_result = {
            "cdr1a_aa": [0,0],
            "cdr2a_aa": [0,0],
            "cdr3a_aa": [0,0],
            "cdr1b_aa": [0,0],
            "cdr2b_aa": [0,0],
            "cdr3b_aa": [0,0],
        }

        if self.train_dataset is not None:
            max_train_sequence = len(self.train_dataset) if max_train_sequence is None else max_train_sequence
            if show_progress:
                pbar = tqdm.tqdm(total=max_train_sequence // n_per_batch)
            for j in range(0, max_train_sequence, n_per_batch):
                epoch_indices = np.array(
                    self.train_dataset[j:j+n_per_batch]['input_ids']
                )

                epoch_attention_mask = np.array(
                    self.train_dataset[j:j+n_per_batch]['attention_mask']
                )

                epoch_indices_mlm, epoch_attention_mask_mlm = self.collator(
                    epoch_indices, epoch_attention_mask
                )
                epoch_indices_mlm = epoch_indices_mlm.to(self.device)
                epoch_attention_mask_mlm = epoch_attention_mask_mlm.to(self.device)

                output = self.model.forward(
                    input_ids = epoch_indices_mlm,
                    attention_mask = epoch_attention_mask_mlm,
                    labels = torch.tensor(
                        epoch_indices, 
                        dtype=torch.int64
                    ).to(self.device)
                )

                predictions = output['output']['logits'].topk(1)[1].squeeze()
                evaluate_mask = epoch_indices_mlm == 4
                ground_truth = torch.tensor(epoch_indices, dtype=torch.int64).to(self.device)

                ground_truth_cdr1a = ground_truth[:,1:9]
                predictions_cdr1a = predictions[:,1:9]
                cdr1a_prediction = (
                    ((ground_truth_cdr1a == predictions_cdr1a) * evaluate_mask[:,1:9].detach()).sum(),
                    evaluate_mask[:,1:9].sum()
                )

                ground_truth_cdr2a = ground_truth[:,9:18]
                predictions_cdr2a = predictions[:,9:18]
                cdr2a_prediction = (
                    ((ground_truth_cdr2a == predictions_cdr2a) * evaluate_mask[:,9:18].detach()).sum(),
                    evaluate_mask[:,9:18].sum()
                )


                ground_truth_cdr3a = ground_truth[:,18:55]
                predictions_cdr3a = predictions[:,18:55]
                cdr3a_prediction = (
                    ((ground_truth_cdr3a == predictions_cdr3a) * evaluate_mask[:,18:55].detach()).sum(),
                    evaluate_mask[:,18:55].sum()
                )

                ground_truth_cdr1b = ground_truth[:,55:64]
                predictions_cdr1b = predictions[:,55:64]
                cdr1b_prediction = (
                    ((ground_truth_cdr1b == predictions_cdr1b) * evaluate_mask[:,55:64].detach()).sum(),
                    evaluate_mask[:,55:64].sum()
                )

                ground_truth_cdr2b = ground_truth[:,64:73]
                predictions_cdr2b = predictions[:,64:73]
                cdr2b_prediction = (
                    ((ground_truth_cdr2b == predictions_cdr2b) * evaluate_mask[:,64:73].detach()).sum(),
                    evaluate_mask[:,64:73].sum()
                )

                ground_truth_cdr3b = ground_truth[:,73:110]
                predictions_cdr3b = predictions[:,73:110]
                cdr3b_prediction = (
                    ((ground_truth_cdr3b == predictions_cdr3b) * evaluate_mask[:,73:110].detach()).sum(),
                    evaluate_mask[:,73:110].sum()
                )
                
                all_train_result['cdr1a_aa'][0] += cdr1a_prediction[0].item()
                all_train_result['cdr1a_aa'][1] += cdr1a_prediction[1].item()
                all_train_result['cdr2a_aa'][0] += cdr2a_prediction[0].item()
                all_train_result['cdr2a_aa'][1] += cdr2a_prediction[1].item()
                all_train_result['cdr3a_aa'][0] += cdr3a_prediction[0].item()
                all_train_result['cdr3a_aa'][1] += cdr3a_prediction[1].item()
                all_train_result['cdr1b_aa'][0] += cdr1b_prediction[0].item()
                all_train_result['cdr1b_aa'][1] += cdr1b_prediction[1].item()
                all_train_result['cdr2b_aa'][0] += cdr2b_prediction[0].item()
                all_train_result['cdr2b_aa'][1] += cdr2b_prediction[1].item()
                all_train_result['cdr3b_aa'][0] += cdr3b_prediction[0].item()
                all_train_result['cdr3b_aa'][1] += cdr3b_prediction[1].item()

                if show_progress:
                    pbar.update(1)
            if show_progress:
                pbar.close()

        ### Test

        all_test_result = {
            "cdr1a_aa": [0,0],
            "cdr2a_aa": [0,0],
            "cdr3a_aa": [0,0],
            "cdr1b_aa": [0,0],
            "cdr2b_aa": [0,0],
            "cdr3b_aa": [0,0],
        }
        if self.test_dataset is not None:
            max_test_sequence = len(self.test_dataset) if max_test_sequence is None else max_test_sequence
            if show_progress:
                pbar = tqdm.tqdm(total=max_test_sequence // n_per_batch)

            for j in range(0, max_test_sequence, n_per_batch):
                epoch_indices = np.array(
                    self.test_dataset[j:j+n_per_batch]['input_ids']
                )

                epoch_attention_mask = np.array(
                    self.test_dataset[j:j+n_per_batch]['attention_mask']
                )

                epoch_indices_mlm, epoch_attention_mask_mlm = self.collator(
                    epoch_indices, epoch_attention_mask
                )
                epoch_indices_mlm = epoch_indices_mlm.to(self.device)
                epoch_attention_mask_mlm = epoch_attention_mask_mlm.to(self.device)

                output = self.model.forward(
                    input_ids = epoch_indices_mlm,
                    attention_mask = epoch_attention_mask_mlm,
                    labels = torch.tensor(
                        epoch_indices, 
                        dtype=torch.int64
                    ).to(self.device)
                )

                predictions = output['output']['logits'].topk(1)[1].squeeze()
                evaluate_mask = epoch_indices_mlm == 4
                ground_truth = torch.tensor(epoch_indices, dtype=torch.int64).to(self.device)

                ground_truth_cdr1a = ground_truth[:,1:9]
                predictions_cdr1a = predictions[:,1:9]
                cdr1a_prediction = (
                    ((ground_truth_cdr1a == predictions_cdr1a) * evaluate_mask[:,1:9].detach()).sum(),
                    evaluate_mask[:,1:9].sum()
                )

                ground_truth_cdr2a = ground_truth[:,9:18]
                predictions_cdr2a = predictions[:,9:18]
                cdr2a_prediction = (
                    ((ground_truth_cdr2a == predictions_cdr2a) * evaluate_mask[:,9:18].detach()).sum(),
                    evaluate_mask[:,9:18].sum()
                )


                ground_truth_cdr3a = ground_truth[:,18:55]
                predictions_cdr3a = predictions[:,18:55]
                cdr3a_prediction = (
                    ((ground_truth_cdr3a == predictions_cdr3a) * evaluate_mask[:,18:55].detach()).sum(),
                    evaluate_mask[:,18:55].sum()
                )

                ground_truth_cdr1b = ground_truth[:,55:64]
                predictions_cdr1b = predictions[:,55:64]
                cdr1b_prediction = (
                    ((ground_truth_cdr1b == predictions_cdr1b) * evaluate_mask[:,55:64].detach()).sum(),
                    evaluate_mask[:,55:64].sum()
                )

                ground_truth_cdr2b = ground_truth[:,64:73]
                predictions_cdr2b = predictions[:,64:73]
                cdr2b_prediction = (
                    ((ground_truth_cdr2b == predictions_cdr2b) * evaluate_mask[:,64:73].detach()).sum(),
                    evaluate_mask[:,64:73].sum()
                )

                ground_truth_cdr3b = ground_truth[:,73:110]
                predictions_cdr3b = predictions[:,73:110]
                cdr3b_prediction = (
                    ((ground_truth_cdr3b == predictions_cdr3b) * evaluate_mask[:,73:110].detach()).sum(),
                    evaluate_mask[:,73:110].sum()
                )

                all_test_result['cdr1a_aa'][0] += cdr1a_prediction[0].item()
                all_test_result['cdr1a_aa'][1] += cdr1a_prediction[1].item()
                all_test_result['cdr2a_aa'][0] += cdr2a_prediction[0].item()
                all_test_result['cdr2a_aa'][1] += cdr2a_prediction[1].item()
                all_test_result['cdr3a_aa'][0] += cdr3a_prediction[0].item()
                all_test_result['cdr3a_aa'][1] += cdr3a_prediction[1].item()
                all_test_result['cdr1b_aa'][0] += cdr1b_prediction[0].item()
                all_test_result['cdr1b_aa'][1] += cdr1b_prediction[1].item()
                all_test_result['cdr2b_aa'][0] += cdr2b_prediction[0].item()
                all_test_result['cdr2b_aa'][1] += cdr2b_prediction[1].item()
                all_test_result['cdr3b_aa'][0] += cdr3b_prediction[0].item()
                all_test_result['cdr3b_aa'][1] += cdr3b_prediction[1].item()

                if show_progress:   
                    pbar.update(1)
            if show_progress:
                pbar.close()
        return all_train_result, all_test_result