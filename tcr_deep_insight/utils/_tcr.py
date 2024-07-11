from dataclasses import dataclass
from typing import Literal
from Bio import SeqIO
from pathlib import Path
import difflib
import numpy as np

from ._amino_acids import _AMINO_ACIDS_INDEX, _AMINO_ACIDS_INDEX_REVERSE

MODULE_PATH = Path(__file__).parent

@dataclass
class MouseTCRAnnotations:
    TRBV_CDR3_REGION_SEQUENCE = {
        'TRBV1': 'CTCSA',
        'TRBV2': 'CASSQ',
        'TRBV3': 'CASSL',
        'TRBV4': 'CASS*',
        'TRBV5': 'CASSQ',
        'TRBV8': 'GASSS',
        'TRBV9': 'CARSL',
        'TRBV10': 'ASSEQ',
        'TRBV12-1': 'CASSL',
        'TRBV12-2': 'CASSL',
        'TRBV13-1': 'CASSD',
        'TRBV13-2': 'CASGD',
        'TRBV13-3': 'CASSD',
        'TRBV14': 'CASSF',
        'TRBV15': 'CASSL',
        'TRBV16': 'CASSL',
        'TRBV17': 'CASSR',
        'TRBV19': 'CASSI',
        'TRBV20': 'CGAR',
        'TRBV21': 'CASSQ',
        'TRBV23': 'CSSSQ',
        'TRBV24': 'CASSL',
        'TRBV26': 'CASSL',
        'TRBV29': 'CASSL',
        'TRBV30': 'CSSR',
        'TRBV31': 'CAWS'
    }

    TRAV_CDR3_REGION_SEQUENCE = {
        'TRAV1': 'CAVR',
        'TRAV2': 'CIVTD',
        'TRAV3-1': 'CAVS',
        'TRAV3-3': 'CAVS',
        'TRAV3-4': 'CAVS',
        'TRAV3D-3': 'CAVS',
        'TRAV3N-3': 'CAVS',
        'TRAV4-2': 'CAAE',
        'TRAV4-3': 'CAAE',
        'TRAV4-4-DV10': 'CAAE',
        'TRAV4D-2': 'CAAE',
        'TRAV4D-3': 'CAAE',
        'TRAV4D-4': 'CAAE',
        'TRAV4N-3': 'CAAE',
        'TRAV4N-4': 'CAAE',
        'TRAV5-1': 'CSAS',
        'TRAV5-2': 'CAES',
        'TRAV5-4': 'CAAS',
        'TRAV5D-2': '*K',
        'TRAV5D-4': 'CAAS',
        'TRAV5N-2': '*K',
        'TRAV5N-4': 'CAAS',
        'TRAV6-1': 'CVLG',
        'TRAV6-2': 'CVLG',
        'TRAV6-3': 'CAMR',
        'TRAV6-4': 'CALV',
        'TRAV6-5': 'CALS',
        'TRAV6-6': 'CALG',
        'TRAV6-7-DV9': 'CALG',
        'TRAV6D-3': 'CAMR',
        'TRAV6D-4': 'CALV',
        'TRAV6D-5': 'CALG',
        'TRAV6D-6': 'CALS',
        'TRAV6D-7': 'CALS',
        'TRAV6N-5': 'CALG',
        'TRAV6N-6': 'CALS',
        'TRAV6N-7': 'CALG',
        'TRAV7-1': 'CAVS',
        'TRAV7-2': 'CAVS',
        'TRAV7-3': 'CAVS',
        'TRAV7-4': 'CAASE',
        'TRAV7-5': 'CAMS',
        'TRAV7-6': 'CAVS',
        'TRAV7D-2': 'CAAS',
        'TRAV7D-3': 'CAVS',
        'TRAV7D-4': 'CAASE',
        'TRAV7D-5': 'CAVS',
        'TRAV7D-6': 'CAAS',
        'TRAV7N-4': 'CAVSE',
        'TRAV7N-5': 'CAVS',
        'TRAV7N-6': 'CAVS',
        'TRAV8-1': 'CATD',
        'TRAV8-2': 'CATD',
        'TRAV8D-1': 'CATD',
        'TRAV8D-2': 'CATD',
        'TRAV8N-2': 'CATD',
        'TRAV9-1': 'CAVS',
        'TRAV9-2': 'CAAS',
        'TRAV9-3': 'CAVS',
        'TRAV9-4': 'CAVS',
        'TRAV9D-1': 'CAAS',
        'TRAV9D-2': 'CAVS',
        'TRAV9D-3': 'CAVS',
        'TRAV9D-4': 'CALS',
        'TRAV9N-2': 'CVLS',
        'TRAV9N-3': 'CAVS',
        'TRAV9N-4': 'CALS',
        'TRAV10': 'CAAS',
        'TRAV10D': 'CAAS',
        'TRAV10N': 'CAAS',
        'TRAV11': 'CVVG',
        'TRAV11D': 'CVVG',
        'TRAV11N': 'CVVG',
        'TRAV12-1': 'CALS',
        'TRAV12-2': 'CALS',
        'TRAV12-3': 'CALS',
        'TRAV12D-1': 'CALS',
        'TRAV12D-2': 'CALS',
        'TRAV12D-3': 'CALS',
        'TRAV12N-1': 'CALS',
        'TRAV12N-2': 'CALS',
        'TRAV12N-3': 'CALS',
        'TRAV13-1': 'CAME',
        'TRAV13-2': 'CAID',
        'TRAV13-3': 'CAME',
        'TRAV13-4-DV7': 'CAME',
        'TRAV13-5': 'CVLS',
        'TRAV13D-1': 'CAME',
        'TRAV13D-2': 'CAID',
        'TRAV13D-3': 'CAME',
        'TRAV13D-4': 'CAME',
        'TRAV13N-1': 'CAME',
        'TRAV13N-2': 'CAID',
        'TRAV13N-3': 'CAME',
        'TRAV13N-4': 'CAME',
        'TRAV14-1': 'CAAS',
        'TRAV14-2': 'CAAS',
        'TRAV14-3': 'CAAS',
        'TRAV14D-1': 'CAAS',
        'TRAV14D-2': 'CAAS',
        'TRAV14D-3-DV8': 'CAAS',
        'TRAV14N-1': 'CAAS',
        'TRAV14N-2': 'CAAS',
        'TRAV14N-3': 'CAAS',
        'TRAV15-1-DV6-1': 'CALWEL',
        'TRAV15-2-DV6-2': 'CALSEL',
        'TRAV15D-1-DV6D-1': 'CALWEL',
        'TRAV15D-2-DV6D-2': 'CALSEL',
        'TRAV15N-1': 'CALWEL',
        'TRAV15N-2': 'CALSEL',
        'TRAV16': 'CAMRE',
        'TRAV16D-DV11': 'CAMRE',
        'TRAV16N': 'CAMRE',
        'TRAV17': 'CALE',
        'TRAV18': 'CAG',
        'TRAV19': 'CAAG',
        'TRAV20': 'CAPG',
        'TRAV21-DV12': 'CILRV'
    }
    TRAJ_CDR3_REGION_SEQUENCE = {
        'TRAJ2': 'NTGGLSGKLTFGEGTQVTVIS',
        'TRAJ3': 'EFSYSSKLIFGAETKLRNPP',
        'TRAJ4': 'LSGSFNKLTFGAGTRLLCAH',
        'TRAJ5': 'GTQVVGQLTFGRGTRLQVYA',
        'TRAJ6': 'TSGGNYKPTFGKGTSLVVHP',
        'TRAJ7': 'DYSNNRLTLGKGTQVVVLP',
        'TRAJ9': 'RNMGYKLTFGTGTSLLVDP',
        'TRAJ11': 'DSGYNKLTFGKGTVLLVSP',
        'TRAJ12': 'GTGGYKVVFGSGTRLLVSP',
        'TRAJ13': 'NSGTYQRFGTGTKLQVVP',
        'TRAJ15': 'YQGGRALIFGTGTTVSVSP',
        'TRAJ16': 'ATSSGQKLVFGQGTILKVYL',
        'TRAJ17': 'TNSAGNKLTFGIGTRVLVRP',
        'TRAJ18': 'DRGSALGRLHFGAGTQLIVIP',
        'TRAJ19': 'IYRGFHKFSSGIESKHNVSP',
        'TRAJ20': 'SGNYKLGVESVTMMSVRA',
        'TRAJ21': 'SNYNVLYFGSGTKLTVEP',
        'TRAJ22': 'SSGSWQLIFGSGTQLTVMP',
        'TRAJ23': ['NYNQGKLIFGQGTKLSIKP', 'PQGKLIF'],
        'TRAJ24': 'ELASLGKLQFGTGTQVVVTP',
        'TRAJ25': 'RTKVSSVFGTWRRLLVKP',
        'TRAJ26': 'NNYAQGLTFGLGTRVSVFP',
        'TRAJ27': 'NTNTGKLTFGDGTVLTVKP',
        'TRAJ28': 'LPGTGSNRLTFGKGTKFSLIP',
        'TRAJ29': 'NSGSRELVLGREARLSMIE',
        'TRAJ30': 'DTNAYKVIFGKGTHLHVLP',
        'TRAJ31': 'NSNNRIFFGDGTQLVVKP',
        'TRAJ32': 'NYGSSGNKLIFGIGTLLSVKP',
        'TRAJ33': ['DSNYQLIWGSGTKLIIKP', 'DGSNYQLIWGSGTKLIIKP'],
        'TRAJ34': ['SSNTDKVVFGTGTRLQVSP', 'SSNTNKVVFGTGTRLQVSP'],
        'TRAJ35': 'QTGFASALTFGSGTKVIPCLP',
        'TRAJ36': 'NIGKKKLVSGTRTRLTIIP',
        'TRAJ37': 'TGNTGKLIFGLGTTLQVQP',
        'TRAJ38': 'NVGDNSKLIWGLGTSLVVNP',
        'TRAJ39': 'NNNAGAKLTFGGGTRLTVRP',
        'TRAJ40': 'VNTGNYKYVFGAGTRLKVIA',
        'TRAJ41': 'VSNTSSMLAEAPHYWSHP',
        'TRAJ42': 'NSGGSNAKLTFGKGTKLSVKS',
        'TRAJ43': 'NNNNAPRFGAGTKLSVKP',
        'TRAJ44': 'VTGSGGKLTLGAGTRLQVNL',
        'TRAJ45': 'NTEGADRLTFGKGTQLIIQP',
        'TRAJ46': 'RRQQCRHAGFGDGDELGVST',
        'TRAJ47': 'HYANKMICGLGTILRVRP',
        'TRAJ48': 'ANYGNEKITFGAGTKLTIKP',
        'TRAJ49': 'NTGYQNFYFGKGTSLTVIP',
        'TRAJ50': 'ASSSFSKLVFGQGTSLSVVP',
        'TRAJ52': 'NTGANTGKLTFGHGTILRVHP',
        'TRAJ53': 'NSGGSNYKLTFGKGTLLTVTP',
        'TRAJ54': 'KRPGFKLVFGQGTGP',
        'TRAJ56': 'ATGGNNKLTFGQGTVLSVIP',
        'TRAJ57': 'NQGGSAKLIFGEGTKLTVSS',
        'TRAJ58': 'QQGTGSKLSFGKGAKLTVSP',
        'TRAJ59': 'LLKREDKATFATGGYEAEED',
        'TRAJ60': 'RSTKDLYFRELSSSSA',
        'TRAJ61': 'VQNEIGIFFFGAMTGRLMKLS'
    }
    TRBJ_CDR3_REGION_SEQUENCE = {
        'TRBJ1-1': 'NTEVFFGKGTRLTVV',
        'TRBJ1-2': 'NSDYTFGSGTRLLVI',
        'TRBJ1-3': 'SGNTLYFGEGSRLIVV',
        'TRBJ1-4': 'SNERLFFGHGTKLSVL',
        'TRBJ1-5': 'NNQAPLFGEGTRLSVL',
        'TRBJ1-6': 'SYNSPLYFAAGTRLTVT',
        'TRBJ1-7': 'PVLDDHGLGKELRYK',
        'TRBJ2-1': 'NYAEQFFGPGTRLTVL',
        'TRBJ2-2': 'NTGQLYFGEGSKLTVL',
        'TRBJ2-3': 'SAETLYFGSGTRLTVL',
        'TRBJ2-4': 'SQNTLYFGAGTRLSVL',
        'TRBJ2-5': 'NQDTQYFGPGTRLLVL',
        'TRBJ2-6': 'ALALTDWQPIEQPMR',
        'TRBJ2-7': 'SYEQYFGPGTRLTVL'
    }
    VJ_GENES = list(TRAV_CDR3_REGION_SEQUENCE.keys()) + list(TRAJ_CDR3_REGION_SEQUENCE.keys()) + list(TRBV_CDR3_REGION_SEQUENCE.keys()) + list(TRBJ_CDR3_REGION_SEQUENCE.keys())

    TRAV_GENES = list(filter(lambda x: x.startswith("TRAV"), VJ_GENES))
    TRAJ_GENES = list(filter(lambda x: x.startswith("TRAJ"), VJ_GENES))
    TRBV_GENES = list(filter(lambda x: x.startswith("TRBV"), VJ_GENES))
    TRBJ_GENES = list(filter(lambda x: x.startswith("TRBJ"), VJ_GENES))

    VJ_GENES2INDEX = dict(zip(VJ_GENES, range(len(VJ_GENES))))
    VJ_GENES2INDEX_REVERSE = {v:k for k,v in VJ_GENES2INDEX.items()}


@dataclass
class HumanTCRAnnotations:
    IMGT_TRAV = [r.id.split('|')[1] for r in SeqIO.parse( MODULE_PATH / "../data/imgt/TRAV.fasta", "fasta")]
    IMGT_TRAJ = [r.id.split('|')[1] for r in SeqIO.parse( MODULE_PATH / "../data/imgt/TRAJ.fasta", "fasta")]
    IMGT_TRBV = [r.id.split('|')[1] for r in SeqIO.parse( MODULE_PATH / "../data/imgt/TRBV.fasta", "fasta")]
    IMGT_TRBD = [r.id.split('|')[1] for r in SeqIO.parse( MODULE_PATH / "../data/imgt/TRBD.fasta", "fasta")]
    IMGT_TRBJ = [r.id.split('|')[1] for r in SeqIO.parse( MODULE_PATH / "../data/imgt/TRBJ.fasta", "fasta")]
    IMGT_VJ_GENES = IMGT_TRAV + IMGT_TRAJ + IMGT_TRBV + IMGT_TRBJ

    TRBV_CDR3_REGION_SEQUENCE = {
        'TRBV1': 'CTSSQ',
        'TRBV2': 'CASSE',
        'TRBV3-1': 'CASSQ',
        'TRBV3-2': 'CASSQ',
        'TRBV4-1': 'CASSQ',
        'TRBV4-2': 'CASSQ',
        'TRBV4-3': 'CASSQ',
        'TRBV5-1': 'CASSL',
        'TRBV5-3': 'CARSL',
        'TRBV5-4': 'CASSL',
        'TRBV5-5': 'CASSL',
        'TRBV5-6': 'CASSL',
        'TRBV5-7': 'CASSL',
        'TRBV5-8': 'CASSL',
        'TRBV6-1': 'CASSE',
        'TRBV6-2': 'CASSY',
        'TRBV6-3': 'CASSY',
        'TRBV6-4': 'CASSD',
        'TRBV6-5': 'CASSY',
        'TRBV6-6': 'CASSY',
        'TRBV6-7': 'CASSY',
        'TRBV6-8': 'CASSY',
        'TRBV6-9': 'CASSY',
        'TRBV7-1': 'CASSS',
        'TRBV7-2': 'CASSL',
        'TRBV7-3': 'CASSL',
        'TRBV7-4': 'CASSL',
        'TRBV7-6': 'CASSL',
        'TRBV7-7': 'CASSL',
        'TRBV7-8': 'CASSL',
        'TRBV7-9': 'CASSL',
        'TRBV9': 'CASSV',
        'TRBV20OR9-2': 'CSAR',
        'TRBV10-1': 'CASSE',
        'TRBV10-2': 'CASSE',
        'TRBV10-3': 'CAISE',
        'TRBV11-1': 'CASSL',
        'TRBV11-2': 'CASSL',
        'TRBV11-3': 'CASSL',
        'TRBV12-1': 'CASSF',
        'TRBV12-2': 'CASRL',
        'TRBV12-3': 'CASSL',
        'TRBV12-4': 'CASSL',
        'TRBV12-5': 'CASGL',
        'TRBV13': 'CASSL',
        'TRBV14': 'CASSQ',
        'TRBV15': 'CATSR',
        'TRBV16': 'CASSQ',
        'TRBV17': 'YSSG',
        'TRBV18': 'CASSP',
        'TRBV19': 'CASSI',
        'TRBV20': 'CSAR',
        'TRBV20-1': 'CSAR',
        'TRBV21-1': 'CASSK',
        'TRBV23-1': 'CASSQ',
        'TRBV24-1': 'CATSDL',
        'TRBV25-1': 'CASSE',
        'TRBV26': 'YASSS',
        'TRBV27': 'CASSL',
        'TRBV28': 'CASSL',
        'TRBV29-1': 'CSVE',
        'TRBV30': 'CAWS'
    }


    TRBJ_CDR3_REGION_SEQUENCE = {
        'TRBJ1-1': 'NTEAFFGQGTRLTVV',
        'TRBJ1-2': 'NYGYTFGSGTRLTVV',
        'TRBJ1-3': 'SGNTIYFGEGSWLTVV',
        'TRBJ1-4': 'TNEKLFFGSGTQLSVL',
        'TRBJ1-5': 'SNQPQHFGDGTRLSIL',
        'TRBJ1-6': 'SYNSPLHFGNGTRLTVT',
        'TRBJ2-1': 'SYNEQFFGPGTRLTVL',
        'TRBJ2-2': 'NTGELFFGEGSRLTVL',
        'TRBJ2-3': 'STDTQYFGPGTRLTVL',
        'TRBJ2-4': 'AKNIQYFGAGTRLSVL',
        'TRBJ2-5': 'QETQYFGPGTRLLVL',
        'TRBJ2-6': 'SGANVLTFGAGSRLTVL',
        'TRBJ2-7': 'SYEQYFGPGTRLTVT'
    }

    TRAV_CDR3_REGION_SEQUENCE = {
    'TRAV1-1': 'CAVR',
    'TRAV1-2': 'CAVR',
    'TRAV2': 'CAVE',
    'TRAV3': 'CAVRD',
    'TRAV4': 'CLVGD',
    'TRAV5': 'CAES',
    'TRAV6': 'CALD',
    'TRAV7': 'CAVD',
    'TRAV8-1': 'CAVN',
    'TRAV8-2': 'CVVS',
    'TRAV8-3': 'CAVG',
    'TRAV8-4': 'CAVS',
    'TRAV8-6': 'CAVS',
    'TRAV8-7': 'CAVG',
    'TRAV9-1': 'CALS',
    'TRAV9-2': 'CALS',
    'TRAV10': 'CVVS',
    'TRAV11': 'CAL',
    'TRAV12-1': 'CVVN',
    'TRAV12-2': 'CAVN',
    'TRAV12-3': 'CAMS',
    'TRAV13-1': 'CAAS',
    'TRAV13-2': 'CAEN',
    'TRAV14DV4': 'CAMRE',
    'TRAV16': 'CALS',
    'TRAV17': 'CATD',
    'TRAV18': 'CALR',
    'TRAV19': 'CALSE',
    'TRAV20': 'CAVQ',
    'TRAV21': 'CAVR',
    'TRAV22': 'CAVE',
    'TRAV23DV6': 'CAAS',
    'TRAV24': 'CAF',
    'TRAV25': 'CAG',
    'TRAV26-1': 'CIVRV',
    'TRAV26-2': 'CILRD',
    'TRAV27': 'CAG',
    'TRAV29DV5': 'CAAS',
    'TRAV30': 'CGTE',
    'TRAV34': 'CGAD',
    'TRAV35': 'CAGQ',
    'TRAV36DV7': 'CAVE',
    'TRAV38-1': 'CAFMK',
    'TRAV38-2DV8': 'CAYRS',
    'TRAV39': 'CAVD',
    'TRAV40': 'CLLG',
    'TRAV41': 'CAVR'
    }


    TRAJ_CDR3_REGION_SEQUENCE = {
        'TRAJ1': 'YESITSQLQFGKGTRVSTSP',
        'TRAJ2': 'NTGGTIDKLTFGKGTHVFIIS',
        'TRAJ3': 'GYSSASKIIFGSGTRLSIRP',
        'TRAJ4': 'FSGGYNKLIFGAGTRLAVHP',
        'TRAJ5': 'DTGRRALTFGSGTRLQVQP',
        'TRAJ6': 'ASGGSYIPTFGRGTSLIVHP',
        'TRAJ7': 'DYGNNRLAFGKGNQVVVIP',
        'TRAJ8': 'NTGFQKLVFGTGTRLLVSP',
        'TRAJ9': 'GNTGGFKTIFGAGTRLFVKA',
        'TRAJ10': 'ILTGGGNKLTFGTGTQLKVEL',
        'TRAJ11': 'NSGYSTLTFGKGTMLLVSP',
        'TRAJ12': 'MDSSYKLIFGSGTRLLVRP',
        'TRAJ13': 'NSGGYQKVTFGIGTKLQVIP',
        'TRAJ14': 'IYSTFIFGSGTRLSVKP',
        'TRAJ15': 'NQAGTALIFGKGTTLSVSS',
        'TRAJ16': 'FSDGQKLLFARGTMLKVDL',
        'TRAJ17': 'IKAAGNKLTFGGGTRVLVKP',
        'TRAJ18': 'DRGSTLGRLYFGRGTQLTVWP',
        'TRAJ19': 'YQRFYNFTFGKGSKHNVTP',
        'TRAJ20': 'SNDYKLSFGAGTTVTVRA',
        'TRAJ21': 'YNFNKFYFGSGTKLNVKP',
        'TRAJ22': 'SSGSARQLTFGSGTQLTVLP',
        'TRAJ23': 'IYNQGGKLIFGQGTELSVKP',
        'TRAJ24': 'TTDSWGKFEFGAGTQVVVTP',
        'TRAJ25': 'EGQGFSFIFGKGTRLLVKP',
        'TRAJ26': 'DNYGQNFVFGPGTRLSVLP',
        'TRAJ27': 'NTNAGKSTFGDGTTLTVKP',
        'TRAJ28': 'YSGAGSYQLTFGKGTKLSVIP',
        'TRAJ29': 'NSGNTPLVFGKGTRLSVIA',
        'TRAJ30': 'NRDDKIIFGKGTRLHILP',
        'TRAJ31': 'NNNARLMFGDGTQLVVKP',
        'TRAJ32': 'NYGGATNKLIFGTGTLLAVQP',
        'TRAJ33': 'DSNYQLIWGAGTKLIIKP',
        'TRAJ34': 'SYNTDKLIFGTGTRLQVFP',
        'TRAJ35': 'IGFGNVLHCGSGTQVIVLP',
        'TRAJ36': 'QTGANNLFFGTGTRLTVIP',
        'TRAJ37': 'GSGNTGKLIFGQGTTLQVKP',
        'TRAJ38': 'NAGNNRKLIWGLGTSLAVNP',
        'TRAJ39': 'NNNAGNMLTFGGGTRLMVKP',
        'TRAJ40': 'TTSGTYKYIFGTGTRLKVLA',
        'TRAJ41': 'NSNSGYALNFGKGTSLLVTP',
        'TRAJ42': 'NYGGSQGNLIFGKGTKLSVKP',
        'TRAJ43': 'NNNDMRFGAGTRLTVKP',
        'TRAJ44': 'NTGTASKLTFGTGTRLQVTL',
        'TRAJ45': 'YSGGGADGLTFGKGTHLIIQP',
        'TRAJ46': 'KKSSGDKLTFGTGTRLAVRP',
        'TRAJ47': 'EYGNKLVFGAGTILRVKS',
        'TRAJ48': 'SNFGNEKLTFGTGTRLTIIP',
        'TRAJ49': 'NTGNQFYFGTGTSLTVIP',
        'TRAJ50': 'KTSYDKVIFGPGTSLSVIP',
        'TRAJ51': 'MRDSYEKLIFGKETLTVKP',
        'TRAJ52': 'NAGGTSYGKLTFGQGTILTVHP',
        'TRAJ53': 'NSGGSNYKLTFGKGTLLTVNP',
        'TRAJ54': 'IQGAQKLVFGQGTRLTINP',
        'TRAJ55': 'KCWCSCWGKGMSTKINP',
        'TRAJ56': 'YTGANSKLTFGKGITLSVRP',
        'TRAJ57': 'TQGGSEKLVFGKGTKLTVNP',
        'TRAJ58': 'ETSGSRLTFGEGTQLTVNP',
        'TRAJ59': 'KEGNRKFTFGMGTQVRVKL',
        'TRAJ60': 'KITMLNFGKGTELIVSL',
        'TRAJ61': 'YRVNRKLTFGANTRGIMKL',
    }

    VJ_GENES = list(TRAV_CDR3_REGION_SEQUENCE.keys()) + list(TRAJ_CDR3_REGION_SEQUENCE.keys()) + list(TRBV_CDR3_REGION_SEQUENCE.keys()) + list(TRBJ_CDR3_REGION_SEQUENCE.keys())

    TRAV_GENES = list(filter(lambda x: x.startswith("TRAV"), VJ_GENES))
    TRAJ_GENES = list(filter(lambda x: x.startswith("TRAJ"), VJ_GENES))
    TRBV_GENES = list(filter(lambda x: x.startswith("TRBV"), VJ_GENES))
    TRBJ_GENES = list(filter(lambda x: x.startswith("TRBJ"), VJ_GENES))

    VJ_GENES2INDEX = dict(zip(VJ_GENES, range(len(VJ_GENES))))
    VJ_GENES2INDEX_REVERSE = {v:k for k,v in VJ_GENES2INDEX.items()}


def get_overlap(s1, s2):
    s = difflib.SequenceMatcher(None, s1, s2)
    pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2))
    return s1[pos_a:pos_a+size]


def getCDR3amr(cdr3a, v, j, species: Literal['human','mouse'] = "human"):
    if species == "human":
        vsequence = HumanTCRAnnotations.TRAV_CDR3_REGION_SEQUENCE[v]
        jsequence = HumanTCRAnnotations.TRAJ_CDR3_REGION_SEQUENCE[j]
    else:
        vsequence = MouseTCRAnnotations.TRAV_CDR3_REGION_SEQUENCE[v]
        jsequence = MouseTCRAnnotations.TRAJ_CDR3_REGION_SEQUENCE[j]
    vmotif = get_overlap(cdr3a, vsequence)
    jmotif = get_overlap(cdr3a, jsequence)
    return cdr3a[cdr3a.index(vmotif)+len(vmotif):cdr3a.index(jmotif)]


def getCDR3bmr(cdr3b, v, j, species: Literal['human','mouse'] = "human"):
    if species == "human":
        vsequence = HumanTCRAnnotations.TRBV_CDR3_REGION_SEQUENCE[v]
        jsequence = HumanTCRAnnotations.TRBJ_CDR3_REGION_SEQUENCE[j]
    else:
        vsequence = MouseTCRAnnotations.TRBV_CDR3_REGION_SEQUENCE[v]
        jsequence = MouseTCRAnnotations.TRBJ_CDR3_REGION_SEQUENCE[j]
    vmotif = get_overlap(cdr3b, vsequence)
    jmotif = get_overlap(cdr3b, jsequence)
    return cdr3b[cdr3b.index(vmotif)+len(vmotif):cdr3b.index(jmotif)]

def getCDR3amr_indices(cdr3a, v, j, species: Literal['human','mouse'] = "human"):
    if species == "human":
        vsequence = HumanTCRAnnotations.TRAV_CDR3_REGION_SEQUENCE[v]
        jsequence = HumanTCRAnnotations.TRAJ_CDR3_REGION_SEQUENCE[j]
    else:
        vsequence = MouseTCRAnnotations.TRAV_CDR3_REGION_SEQUENCE[v]
        jsequence = MouseTCRAnnotations.TRAJ_CDR3_REGION_SEQUENCE[j]
    vmotif = get_overlap(cdr3a, vsequence)
    jmotif = get_overlap(cdr3a, jsequence)
    return cdr3a.index(vmotif)+len(vmotif),cdr3a.index(jmotif)


def getCDR3bmr_indices(cdr3b, v, j, species: Literal['human','mouse'] = "human"):
    if species == "human":
        vsequence = HumanTCRAnnotations.TRBV_CDR3_REGION_SEQUENCE[v]
        jsequence = HumanTCRAnnotations.TRBJ_CDR3_REGION_SEQUENCE[j]
    else:
        vsequence = MouseTCRAnnotations.TRBV_CDR3_REGION_SEQUENCE[v]
        jsequence = MouseTCRAnnotations.TRBJ_CDR3_REGION_SEQUENCE[j]
    vmotif = get_overlap(cdr3b, vsequence)
    jmotif = get_overlap(cdr3b, jsequence)
    return cdr3b.index(vmotif)+len(vmotif),cdr3b.index(jmotif)

def encode_mr_mask(cdr3a, cdr3b, va_gene, ja_gene, vb_gene, jb_gene, length=48, species="human"):
    a = np.zeros(length)
    i,j = getCDR3amr_indices(cdr3a, va_gene, ja_gene, species=species)
    a[i+2:j+2] = 1
    b = np.zeros(length)
    i,j = getCDR3bmr_indices(cdr3b, vb_gene, jb_gene, species=species)
    b[i+2:j+2] = 1
    return np.hstack([a,b])

def encode_tcr(
    annotations,
    trav,
    cdr3a,
    traj,
    trbv,
    cdr3b,
    trbj,
    cdr3a_max_length=45,
    cdr3b_max_length=45,
):
    l = len(_AMINO_ACIDS_INDEX)
    tcr_num_array = [annotations.VJ_GENES2INDEX(trav)+l+1] + \
        list(map(lambda x: _AMINO_ACIDS_INDEX[x]+1, cdr3a)) + [0] * (cdr3a_max_length - len(cdr3a)) + \
        [annotations.VJ_GENES2INDEX(traj)+l+1] + \
        [annotations.VJ_GENES2INDEX(trbv)+l+1] + \
        list(map(lambda x: _AMINO_ACIDS_INDEX[x]+1, cdr3b)) + [0] * (cdr3b_max_length - len(cdr3b)) + \
        [annotations.VJ_GENES2INDEX(trbj)+l+1]
    