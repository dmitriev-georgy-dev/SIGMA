"""Локальная база данных хэшей и последовательностей белков для SIGMA v1.1 (ИСПРАВЛЕННАЯ)"""
import hashlib

def _generate_hash(sequence: str) -> str:
    """Вспомогательная функция для генерации MD5-хэша последовательности"""
    return hashlib.md5(sequence.encode()).hexdigest()

# Полный реестр доступных белков с их названиями и последовательностями ДНК
PROTEIN_REGISTRY = {
    "Человеческий проинсулин (1ZEG)": {
        "dna": "tttgtgaaccaacacctgtgcggctcacacctggtggaagctctctacctagtgtgcggggaacgaggcttcttctacacacccaagacccgccgggaggcagaggacctgcaggtggggcaggtggagctgggcgggggccctggtgcaggcagcctgcagcccttggccctggaggggtccctgcagaagcgtggcattgtggaacaatgctgtaccagcatctgctccctctaccagctggagaactactgcaactag",
        "file": "pdb1zeg.ent"
    },
    "Лизоцим куриного яйца (1LYZ)": {
        "dna": "aaagtctttggacgatgtgagctggcagcggctatgaagagtcacggacttgataactaccgcggatacagcctgggaaactgggtgtgtgctgcaaaattcgagagtaacttcaacacccaggctacaaaccgtaacaccgatggcagtaccgactacggcatcctacagatcaacagccgctggtggtgcaacgatggccgaaccccaggatcgaggaacctgtgcaacatcccgtgctcagccctgctgagctcagacataacagcgagcgtgaactgcgcgaagaagatcgtcagcgatggaaacggcatgaacgcgtgggtcgcctggcgcaaccgctgcaagggtaccgacgtccaggcgtggatcagaggctgccggctgtaa",
        "file": "pdb1lyz.ent"
    },
    "Миоглобин кашалота (1MBO)": {
        "dna": "atgggtttgacggatggcgaatggcagctggtgctgaacgtgtggggcaaggtggaggccgacattccgggccacggccaggaggtgctgatccgcctgttcaagggccacccggagaccctggagaagttcgacaagttcaagcacctgaagagcgaggacgagatgaaggccagcgaggacctgaagaagcacggctgcaccgtgctgaccgccctgggcaccatcctgaagaagaagggccagcacgccgccgagatccagccgctggcccagagccacgccaccaagcacaagatcccggtgaagtacctggagttcatcagcgaggccatcatccaggtgctgcagagcaagcacccgggcgacttcggcgccgacgcccagggcgccatgaacaaggccctggagctgttccgcaaggacatcgccgccaagtacaaggagctgggctaccagggctag",
        "file": "pdb1mbo.ent"
    },
    "Человеческий гемоглобин — Альфа-цепь (1A3N)": {
        "dna": "gtgctgtctcccgccgacaagaccaacgtcaaggccgcctggggtaaggtcggcgcgcacgctggcgaggtgggcgccgaggccctggagaggatgttcctgtccttccccaccaccaagacctacttcccgcacttcgacctgagccacggctctgcccaggttaagggccacggcaagaaggtggccgacgcgctgaccaacgccgtggcgcacgtggacgacatgcccaacgcgctgtccgccctgagcgacctgcacgcgcacaagcttcgggtggacccggtcaacttcaagctcctaagccactgcctgctggtgaccctggccgcccacctccccgccgagttcacccctgcggtccacgcctccctggacaagttcctggcttctgtgagcaccgtgctgacctccaaataccgt",
        "file": "pdb1a3n.ent"
    },
    "Зеленый флуоресцентный белок GFP (1EMA)": {
        "dna": "atggtgagcaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaactacaacagccacaacgtctatatcatggccgacaagcagaagaacggcatcaaggtgaacttcaagatccgccacaacatcgaggacggcagcgtgcagctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccgccctgagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactctcggcatggacgagctgtacaag",
        "file": "pdb1ema.ent"
    },
    "Цитохром C сердца лошади (1HRC)": {
        "dna": "ggtgatgttgagaagggcaagaagatctttgtccagaagtgtgcccagtgccacactgtggaaaagggaggcaagcacaagactgggccaaacctccatggtctgtttgggcggaagacagggcaggcccctggatttacttacacagatgccaacaagaacaagggcatcacctggaaagaggagacactgatggagtatttggagaatcccaagaagtacatccctggcacaaaaatgatctttgctggcattaagaagaagacggagcgggaggacctgatcgcttacctcaagaaggctaccaatgagtga",
        "file": "pdb1hrc.ent"
    }
}

# Словарь: хэш аминокислотной последовательности -> имя файла в папке structures
LOCAL_PROTEIN_DB = {
    _generate_hash("FVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"): "pdb1zeg.ent",
    _generate_hash("KVFGRCELAAAMKSHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL"): "pdb1lyz.ent",
    _generate_hash("MGLTDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGCTVLTALGTILKKKGQHAAEIQPLAQSHATKHKIPVKYLEFISEAIIQVLQSKHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG"): "pdb1mbo.ent",
    _generate_hash("VLSPADKTNVKAAWGKVGAHAGEVGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR"): "pdb1a3n.ent",
    _generate_hash("MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"): "pdb1ema.ent",
    _generate_hash("GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"): "pdb1hrc.ent"
}
