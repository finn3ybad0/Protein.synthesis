#Name: Finney Bado
#Date: 2024-03-08
#Description: This program synthetise a protein amino-acid sequence given a DNA nucleotid sequence
#The result is displayed as grid


#######################################################################################################################
#                                        DONNÉES
#######################################################################################################################

adn="TCGACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGACG\
GCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGAGTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGAATGCCAGCCAGC\
CAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGAACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATATTCAGCCAGCCAGCCAGCGA\
ACTCGACACTCTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCATCCCAGCGATACCC\
AGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAG\
CCAGCGAACTCGTCTGCGTTCGACAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGC\
GATTGCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATTGCCAGCCAGCCAGCCAGCGAACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGTATGCCAGCC\
AGCATCCCAGCGA"

codons_aa = {
    "UUU": "Phénylalanine",
    "UUC": "Phénylalanine",
    "UUA": "Leucine",
    "UUG": "Leucine",
    "CUU": "Leucine",
    "CUC": "Leucine",
    "CUA": "Leucine",
    "CUG": "Leucine",
    "AUU": "Isoleucine",
    "AUC": "Isoleucine",
    "AUA": "Isoleucine",
    "AUG": "Méthionine (Start)",
    "GUU": "Valine",
    "GUC": "Valine",
    "GUA": "Valine",
    "GUG": "Valine",
    "UCU": "Sérine",
    "UCC": "Sérine",
    "UCA": "Sérine",
    "UCG": "Sérine",
    "CCU": "Proline",
    "CCC": "Proline",
    "CCA": "Proline",
    "CCG": "Proline",
    "ACU": "Thrénine",
    "ACC": "Thrénine",
    "ACA": "Thrénine",
    "ACG": "Thrénine",
    "GCU": "Alanine",
    "GCC": "Alanine",
    "GCA": "Alanine",
    "GCG": "Alanine",
    "UAU": "Tyrosine",
    "UAC": "Tyrosine",
    "UAA": "Stop",
    "UAG": "Stop",
    "CAU": "Histidine",
    "CAC": "Histidine",
    "CAA": "Glutamine",
    "CAG": "Glutamine",
    "AAU": "Asparagine",
    "AAC": "Asparagine",
    "AAA": "Lysine",
    "AAG": "Lysine",
    "GAU": "Aspartate",
    "GAC": "Aspartate",
    "GAA": "Glutamate",
    "GAG": "Glutamate",
    "UGU": "Cystéine",
    "UGC": "Cystéine",
    "UGA": "Stop",
    "UGG": "Tryptophane",
    "CGU": "Arginine",
    "CGC": "Arginine",
    "CGA": "Arginine",
    "CGG": "Arginine",
    "AGU": "Sérine",
    "AGC": "Sérine",
    "AGA": "Arginine",
    "AGG": "Arginine",
    "GGU": "Glycine",
    "GGC": "Glycine",
    "GGA": "Glycine",
    "GGG": "Glycine"}

lettreAa = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "UGU": "C",
    "UGC": "C",
    "UGA": "*",
    "UGG": "W",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G"
}


#######################################################################################################################
#                                        SOUS-FONCTIONS ET TESTS UNITAIRES
#######################################################################################################################



def antisens(brinAdn):
    if brinAdn == "":
        print('invalid input : no sequence')
    else :
        complement={"A":"T","G":"C","T":"A","C":"G"}
        brinCompl=""
        for i in brinAdn:
            if i not in "ATGC" :
                print("invalid input : sequence contains invalid nucleotid")
                break
            else:
                brinCompl= complement[i] + brinCompl
        return brinCompl


def testAntisens(brinAdn):
    pass  #take in account invalid nucleotid and also empty input

def trouveDebut(brinAdn):
    position=[]
    index = brinAdn.find("TAC")
    while index != -1 :
        position.append(index)
        index = brinAdn.find("TAC", index + 3)
    return position


def testTrouveDebut(brinAdn):
    pass

def trouveFin(brinAdn):
    position = []
    indexATT = brinAdn.find("ATT")
    indexACT = brinAdn.find("ACT")
    indexATC = brinAdn.find("ATC")
    while indexATT!=-1 or indexACT!= -1 or indexATC != -1:
        if indexATT!=-1:
            position.append(indexATT)
            indexATT = brinAdn.find("ATT", indexATT + 3)
        if indexACT!=-1:
            position.append(indexACT)
            indexACT = brinAdn.find("ACT",indexACT + 3)
        if indexATC!=-1:
            position.append(indexATC)
            indexATC = brinAdn.find("ATC", indexATC + 3)
    position.sort()
    return position

#
# def testTrouveFin(brinAdn):
#     assert
#     assert
#     assert
#     assert
#
def trouveGene(debut, fin):
    pointerDebut=0
    pointerFin=0
    positionGene=[]
    while pointerFin < len(fin) and pointerDebut < len(debut):
        if debut[pointerDebut] < fin[pointerFin]:
            if (fin[pointerFin] - debut[pointerDebut] )% 3 == 0 :
                positionGene.append((debut[pointerDebut],fin[pointerFin]))
                pointerDebut += 1
        pointerFin +=1
    return positionGene
#
# def testTrouveGene(debut,fin):
#     assert
#     assert
#     assert
#     assert
#
def transcrire(brinAdn):
    complementArn={"A":"U","G":"C","T":"A","C":"G"}
    brinArn=""
    for i in brinAdn :
        brinArn+= complementArn[i]
    return brinArn

#     pass
#
# def testTranscrire(brinAdn):
#     assert
#     assert
#     assert
#     assert
#
def traduire(brinArn):
    pointer=0
    protein=""
    while pointer<len(brinArn)-3:
            print(pointer)
            codon= brinArn[pointer] + brinArn[pointer+1] + brinArn[pointer+2]
            if pointer == len(brinArn)-6:
                protein += codons_aa[codon]
            else :
                protein += codons_aa[codon] + " - "
            pointer+=3
    return protein

def testTraduire(brinArn):
    pass
def carre(longueur, nombre):
    pass

def testCarre(longueur,nombre):
    pass








#######################################################################################################################
#                                        FONCTION PRINCIPALE ET TEST UNITAIRES
#######################################################################################################################


def adnToProtein(brinAdn):
    adn=[brinAdn,antisens(brinAdn)]
    print(adn)
    genes=[]
    for i in adn:
        genes.append(trouveGene(trouveDebut(i),trouveFin(i)))
    # for i in genes:
    return genes

