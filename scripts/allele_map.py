def get_proper_name(allele: str):
    if allele.lower() == 'y':
        return 'c/t'
    elif allele.lower() == 'r':
        return 'a/g'
