# Hugginface Transformers
from transformers import BertConfig

from ...utils._compat import Literal

BERT_TINY = {
    "hidden_size": 128,
    "intermediate_size": 512,
    "num_hidden_layers": 4, 
    "num_attention_heads": 4
}

BERT_SMALL = {
    "hidden_size": 192,
    "intermediate_size": 768,
    "num_hidden_layers": 6, 
    "num_attention_heads": 6
}


def get_config(
    extra_vocab_size=183,
    **kwargs
):
    bert_base = {
        "attention_probs_dropout_prob": 0.1,
        "hidden_act": "gelu",
        "hidden_dropout_prob": 0.1,
        "hidden_size": 768,
        "initializer_range": 0.02,
        "intermediate_size": 1536,
        "layer_norm_eps": 1e-12,
        "max_position_embeddings": 516,
        "model_type": "bert",
        "num_attention_heads": 12,
        "num_hidden_layers": 12,
        "pad_token_id": 0,
        "position_embedding_type": "absolute",
        "type_vocab_size": 2,
        "use_cache": 1,
        "vocab_size": 27+extra_vocab_size # Amino Acids + V-J genes
    }
    bert_base.update(kwargs)
    return bert_base


def get_human_config(
    bert_type: Literal['tiny','small','base']="small", 
    vocab_size: int = None,
    alibi_starting_size: int = 512,
) -> BertConfig:
    """
    Get the configuration for the human TCR BERT

    :param bert_type: The size of the BERT model. Must be one of 'tiny', 'small', or 'base'
    :param vocab_size: The size of the vocabulary
    :param alibi_starting_size: The size of the input sequence

    :return: The configuration for the human TCR BERT
    """
    if bert_type == "tiny":
        config = get_config(**BERT_TINY)
    elif bert_type == "small":
        config = get_config(**BERT_SMALL)
    elif bert_type == "base":
        config = get_config()
    else:
        raise ValueError("bert_type must be one of 'tiny', 'small', or 'base'")
    
    if vocab_size is not None:
        config["vocab_size"] = vocab_size
    config['alibi_starting_size'] = alibi_starting_size
    config['attention_probs_dropout_prob'] = 0.0
    return BertConfig.from_dict(config)

def get_mouse_config(
    bert_type: Literal['tiny','small','base']="small", 
    vocab_size: int = None,
    alibi_starting_size: int = 512,
) -> BertConfig:
    """
    Get the configuration for the mouse TCR BERT

    :param bert_type: The size of the BERT model. Must be one of 'tiny', 'small', or 'base'
    :param vocab_size: The size of the vocabulary
    :param alibi_starting_size: The size of the input sequence

    :return: The configuration for the mouse TCR BERT
    """
    if bert_type == "tiny":
        return BertConfig.from_dict(get_config(extra_vocab_size=213, **BERT_TINY, ))
    elif bert_type == "small":
        return BertConfig.from_dict(get_config(extra_vocab_size=213, **BERT_SMALL))
    elif bert_type == "base":
        return BertConfig.from_dict(get_config(extra_vocab_size=213))
    else:
        raise ValueError("bert_type must be one of 'tiny', 'small', or 'base'")