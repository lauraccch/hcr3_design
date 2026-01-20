def amp(ampl):
    if ampl == "B1":
        upspc="aa"; dnspc="ta"; up="GAGGAGGGCAGCAAACGG"; dn="GAAGAGTCTTCCTTTACG"
    elif ampl == "B2":
        upspc="aa"; dnspc="aa"; up="CCTCGTAAATCCTCATCA"; dn="ATCATCCAGTAAACCGCC"
    elif ampl == "B3":
        upspc="tt"; dnspc="tt"; up="GTCCCTGCCTCTATATCT"; dn="CCACTCAACTTTAACCCG"
    elif ampl == "B4":
        upspc="aa"; dnspc="at"; up="CCTCAACCTACCTCCAAC"; dn="TCTCACCATATTCGCTTC"
    elif ampl == "B5":
        upspc="aa"; dnspc="aa"; up="CTCACTCCCAATCTCTAT"; dn="CTACCCTACAAATCCAAT"
    elif ampl == "B7":
        upspc="ww"; dnspc="ww"; up="CTTCAACCTCCACCTACC"; dn="TCCAATCCCTACCCTCAC"
    elif ampl == "B9":
        upspc="ww"; dnspc="ww"; up="CACGTATCTACTCCACTC"; dn="TCAGCACACTCCCAACCC"
    elif ampl == "B10":
        upspc="ww"; dnspc="ww"; up="CCTCAAGATACTCCTCTA"; dn="CCTACTCGACTACCCTAG"
    elif ampl == "B11":
        upspc="ww"; dnspc="ww"; up="CGCTTAGATATCACTCCT"; dn="ACGTCGACCACACTCATC"
    elif ampl == "B13":
        upspc="ww"; dnspc="ww"; up="AGGTAACGCCTTCCTGCT"; dn="TTATGCTCAACATACAAC"
    elif ampl == "B14":
        upspc="ww"; dnspc="ww"; up="AATGTCAATAGCGAGCGA"; dn="CCCTATATTTCTGCACAG"
    elif ampl == "B15":
        upspc="ww"; dnspc="ww"; up="CAGATTAACACACCACAA"; dn="GGTATCTCGAACACTCTC"
    elif ampl == "B17":
        upspc="ww"; dnspc="ww"; up="CGATTGTTTGTTGTGGAC"; dn="GCATGCTAATCGGATGAG"
    else:
        raise ValueError("Unknown amplifier type")
    return [upspc, dnspc, up, dn]

BOLD = '\033[1m'
END = '\033[0m'
tab = "\t"