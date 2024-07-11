#!/usr/bin/env python3
import pandas as pd 

huardb_annotation_low_cmap = {
    'CD4': '#1c83c5ff',
    'CD8': '#ffca39ff',
    'CD40LG': '#5bc8d9ff',
    'Cycling': '#a7a7a7ff',
    'MAIT': '#2a9d8fff',
    'Naive CD4': '#3c3354ff',
    'Naive CD8': '#a9d55dff',
    'Treg': '#6a4d93ff',
    'Undefined': '#f7f7f7ff',
    'Ambiguous': '#f7f7f7ff',
    'Unknown': '#f7f7f7ff'
}

huardb_annotation_high_cmap = {
    'CD8+ Tex': '#7131D5',
    'CD8+ Tm': '#AA4053',
    'CD4+ Treg': '#AA40FF',
    'CD4+ Tn': '#FFDB8D',
    'CD4+ Tm': '#FF858D',
    'CD8+ Trm': '#0FAF8D',
    'Unpredicted': '#A7A7A7',
    'CD8+ Teff': '#997276',
    'CD8+ Tpex': '#F59495',
    'CD8+ Temra': '#F2C96D',
    'CD8+ Tn': '#CEBF8F',
    'CD8+ Tcm': '#ffbb78',
    'CD8+ Early Tcm/Tem': '#ff7f0e',
    'CD8+ GZMK+ Tem': '#d62728',
    'CD8+ CREM+ Tm': '#aa40fc',
    'CD8+ KLRG1+ Temra': '#8c564b',
    'CD8+ KLRG1- Temra': '#e377c2',
    'CD8+ IFITM3+KLRG1+ Temra': '#b5bd61',
    'CD8+ MAIT': '#7FCA00',
    'CD8+ ILTCK': '#aec7e8',
    'CD8+ ITGAE+ Trm': '#279e68',
    'CD8+ ITGB2+ Trm': '#98df8a',
    'CD8+ SELL+ progenitor Tex': '#ff9896',
    'CD8+ GZMK+ Tex': '#c5b0d5',
    'CD8+ CXCR6+ Tex': '#c49c94',
    'CD8+ Cycling T': '#f7b6d2',
    'CD8+ CD45hi IL32+ Tem': '#00A1A1',
    'CD8+ CD45hi Temra': '#00CACA',
    'CD8+ CD45hi Tn': '#FFD13F',
    'CD8+ CD45hi Trm': '#004F4F',
    'CD8+ IL32+ MIF+ T': '#AA7840',
    'CD8+ IL32+ GZMB+ Teff': '#CF89CF',
    'CD8+ IL32+ TIGIT+ Tex': '#991481',
    'CD4+ Tem': '#D20F8C',
    'CD4+ FOS+ Tn': '#DBDB8D',
    'CD4+ Tcm': '#BF4151',
    'CD4+ KLF2+ Tn': '#CD912D',
    'CD4+ GIMAP4+ Tn': '#D3BD2A',
    'CD4+ ICOS+ Treg': '#1ABCDC',
    'CD4+ FOS+ Tm': '#B445BA',
    'CD4+ Naive-like Treg': '#008FB5',
    'CD4+ AREG+ Tm': '#4653B7',
    'CD4+ IFNG+ T': '#B445BA',
    'CD4+ CCL5+ Tcm': '#E83A4F',
    'CD4+ Tfh': '#A89CFF',
    'CD4+ Cytotoxic T': '#773115',
    'CD4+ Cycling T': '#F7B6D1',
    'CD4+ CCR7+ Tm': '#FA8115',
    'CD4+ CXCR5+ Treg': '#4D3E98',
    'CD4': '#1c83c5',
    'CD8': '#ffca39',
    'CD40LG': '#5bc8d9',
    'Cycling': '#a7a7a7',
    'MAIT': '#2a9d8f',
    'Naive CD4': '#3c3354',
    'Naive CD8': '#a9d55d',
    'Treg': '#6a4d93',
    'Undefined': '#f7f7f7',
    'Ambiguous': '#f7f7f7',
    'Unknown': '#f7f7f7'
}


godsnot_102 = [
    '#FFFF00', '#1CE6FF', '#FF34FF', '#FF4A46', '#008941', '#006FA6',
    '#A30059', '#FFDBE5', '#7A4900', '#0000A6', '#63FFAC', '#B79762',
    '#004D43', '#8FB0FF', '#997D87', '#5A0007', '#809693', '#6A3A4C',
    '#1B4400', '#4FC601', '#3B5DFF', '#4A3B53', '#FF2F80', '#61615A',
    '#BA0900', '#6B7900', '#00C2A0', '#FFAA92', '#FF90C9', '#B903AA',
    '#D16100', '#DDEFFF', '#000035', '#7B4F4B', '#A1C299', '#300018',
    '#0AA6D8', '#013349', '#00846F', '#372101', '#FFB500', '#C2FFED',
    '#A079BF', '#CC0744', '#C0B9B2', '#C2FF99', '#001E09', '#00489C',
    '#6F0062', '#0CBD66', '#EEC3FF', '#456D75', '#B77B68', '#7A87A1',
    '#788D66', '#885578', '#FAD09F', '#FF8A9A', '#D157A0', '#BEC459',
    '#456648', '#0086ED', '#886F4C', '#34362D', '#B4A8BD', '#00A6AA',
    '#452C2C', '#636375', '#A3C8C9', '#FF913F', '#938A81', '#575329',
    '#00FECF', '#B05B6F', '#8CD0FF', '#3B9700', '#04F757', '#C8A1A1',
    '#1E6E00', '#7900D7', '#A77500', '#6367A9', '#A05837', '#6B002C',
    '#772600', '#D790FF', '#9B9700', '#549E79', '#FFF69F', '#201625',
    '#72418F', '#BC23FF', '#99ADC0', '#3A2465', '#922329', '#5B4534',
    '#FDE8DC', '#404E55', '#0089A3', '#CB7E98', '#A4E804', '#324E72'
]

default_10 = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
    '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
]

default_20 = [
    '#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
    '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
    '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5',
    '#ad494a', '#8c6d31'
]

default_28 = [
    '#023fa5', '#7d87b9', '#bec1d4', '#d6bcc0', '#bb7784', '#8e063b',
    '#4a6fe3', '#8595e1', '#b5bbe3', '#e6afb9', '#e07b91', '#d33f6a',
    '#11c638', '#8dd593', '#c6dec7', '#ead3c6', '#f0b98d', '#ef9708',
    '#0fcfc0', '#9cded6', '#d5eae7', '#f3e1eb', '#f6c4e1', '#f79cd4',
    '#7f7f7f', '#c7c7c7', '#1CE6FF', '#336600'
]


zheng_2020_annotation_cmap_cd8 = {'CD8.c01.Tn.MAL': '#96C3D8',
 'CD8.c02.Tm.IL7R': '#5D9BBE',
 'CD8.c03.Tm.RPS12': '#F5B375',
 'CD8.c04.Tm.CD52': '#C0937E',
 'CD8.c05.Tem.CXCR5': '#67A59B',
 'CD8.c06.Tem.GZMK': '#A4D38E',
 'CD8.c07.Temra.CX3CR1': '#4A9D47',
 'CD8.c08.Tk.TYROBP': '#F19294',
 'CD8.c09.Tk.KIR2DL4': '#E45A5F',
 'CD8.c10.Trm.ZNF683': '#3477A9',
 'CD8.c11.Tex.PDCD1': '#BDA7CB',
 'CD8.c12.Tex.CXCL13': '#684797',
 'CD8.c13.Tex.myl12a': '#9983B7',
 'CD8.c14.Tex.TCF7': '#CD9A99',
 'CD8.c15.ISG.IFIT1': '#DD4B52',
 'CD8.c16.MAIT.SLC4A10': '#DA8F6F',
 'CD8.c17.Tm.NME1': '#F58135'}

zheng_2020_annotation_cmap_cd4 = {'CD4.c01.Tn.TCF7': '#78AECB',
 'CD4.c02.Tn.PASK': '#639FB0',
 'CD4.c03.Tn.ADSL': '#98C7A5',
 'CD4.c04.Tn.il7r': '#83C180',
 'CD4.c05.Tm.TNF': '#B2A4A5',
 'CD4.c06.Tm.ANXA1': '#EC8D63',
 'CD4.c07.Tm.ANXA2': '#CFC397',
 'CD4.c08.Tm.CREM': '#F6B279',
 'CD4.c09.Tm.CCL5': '#6197B4',
 'CD4.c10.Tm.CAPG': '#CEA168',
 'CD4.c11.Tm.GZMA': '#A0A783',
 'CD4.c12.Tem.GZMK': '#9ACC90',
 'CD4.c13.Temra.CX3CR1': '#6A9A52',
 'CD4.c14.Th17.SLC4A10': '#E97679',
 'CD4.c15.Th17.IL23R': '#DE4247',
 'CD4.c16.Tfh.CXCR5': '#A38CBD',
 'CD4.c17.TfhTh1.CXCL13': '#795FA3',
 'CD4.c18.Treg.RTKN2': '#E0C880',
 'CD4.c19.Treg.S1PR1': '#C28B65',
 'CD4.c20.Treg.TNFRSF9': '#A65A34',
 'CD4.c21.Treg.OAS1': '#DE4B3F',
 'CD4.c22.ISG.IFIT1': '#DD9E82',
 'CD4.c23.Mix.NME1': '#E78B75',
 'CD4.c24.Mix.NME2': '#F7A96C',
 'undefined': '#FFFFFF'}

zheng_2020_annotation_cmap = zheng_2020_annotation_cmap_cd8.copy()
zheng_2020_annotation_cmap.update(zheng_2020_annotation_cmap_cd4)

_chu_2023_annotation_string = """
CD8-3	CD8_c3_Tn	#E9ADC2
CD8-13	CD8_c13_Tn_TCF7	#AACC65
CD8-0	CD8_c0_Teff	#00AFCA
CD8-2	CD8_c2_Teff	#BBB7CB
CD8-8	CD8_c8_Teff_KLRG1	#E1A276
CD8-10	CD8_c10_Teff_CD244	#A5A2B3
CD8-11	CD8_c11_Teff_SEMA4A	#A3AFA9
CD8-6	CD8_c6_Tcm	#DD7A80
CD8-12	CD8_c12_Trm	#A4BD83
CD8-7	CD8_c7_Tpex	#EB9B7F
CD8-1	CD8_c1_Tex	#76BCD8
CD8-4	CD8_c4_Tstr	#E27C97
CD8-5	CD8_c5_Tisg	#DF6C87
CD8-9	CD8_c9_Tsen	#CCA891
CD4-2	CD4_c2_Tn	#E0C8D9
CD4-6	CD4_c6_Tn_FHIT	#F0A683
CD4-7	CD4_c7_Tn_TCEA3	#E5AE7C
CD4-9	CD4_c9_Tn_TCF7_SLC40A1	#A6AEBE
CD4-10	CD4_c10_Tn_LEF1_ANKRD55	#B3C28B
CD4-0	CD4_c0_Tcm	#4CBBD2
CD4-5	CD4_c5_CTL	#E9949E
CD4-1	CD4_c1_Treg	#9FC6DB
CD4-3	CD4_c3_TFH	#EFB3CC
CD4-8	CD4_c8_Th17	#C4ADA6
CD4-4	CD4_c4_Tstr	#EF9AB9
CD4-11	CD4_c11_Tisg	#C4D960
"""

_chu_2023_annotation = _chu_2023_annotation_string.split("\n")[1:-1]
_chu_2023_annotation = list(map(lambda x: x.split('\t'), _chu_2023_annotation))
_chu_2023_annotation = pd.DataFrame(_chu_2023_annotation)
chu_2023_annotation_cmap = dict(zip(_chu_2023_annotation.iloc[:,1], _chu_2023_annotation.iloc[:,2]))
