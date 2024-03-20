# Name: Finney Bado & Johnny Kana
# Date: 2024-03-08
# Description: This program synthesize a protein amino-acid sequence given a DNA nucleotide sequence
# The result is displayed as grid


#######################################################################################################################
#                                        DONNÉES
#######################################################################################################################

import turtle


adn = "TCGACTGCGATCGACAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGAAGCCAGCCAGCCGATACCCAGCCAGCCAGCCAGCGACG\
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
    "GGG": "Glycine",
}

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
    "GGG": "G",
}


#######################################################################################################################
#                                        SOUS-FONCTIONS ET TESTS UNITAIRES
#######################################################################################################################


def antisens(brinAdn: str) -> str:

    # define a dict that maps every nucleotide to it's complement
    complement = {
        "A": "T",
        "G": "C",
        "T": "A",
        "C": "G"
    }
    # Initialize an empty string to store the complement of the DNA strand
    brinCompl = ""

    # concatenate every nucleotide to the beginning of the complement string so that they can be read left to right
    for i in brinAdn:
        brinCompl = complement[i] + brinCompl

    return brinCompl


def testAntisens(brinAdn):

    assert antisens("AAAAAAAAAAAAAAA") == "TTTTTTTTTTTTTTT"
    assert antisens("TTTTTTTTTTTTTTT") == "AAAAAAAAAAAAAAA"
    assert antisens("ATCGATCGATCGATC") == "CGATCGATCGATCGA"


def trouveDebut(brinAdn: str) -> list:

    # Initialize an empty list to store the positions of gene strand beginnings
    position = []
    # Find the position of the first occurrence of "TAC" sequence in the DNA strand (if it exists )
    index = brinAdn.find("TAC")

    # Find the positions of the next occurrences of "TAC" and add them to the list (if they exist)
    while index != -1:
        position.append(index)
        index = brinAdn.find("TAC", index + 3)

    return position


def testTrouveDebut(brinAdn):

    assert trouveDebut("ACGTCGTAGCGATCGATCGTACGTCGTAGCTAG") == []
    assert trouveDebut("TACGTCGATCGTACGTAGCGATCGATCGTAGCT") == [0, 11]
    assert trouveDebut("ACGTCGATCGTATAGCTATAGCTAGCTAGGTAC") == [33]


def trouveFin(brinAdn: str) -> list:
    # Initialize
    position = []
    indexATT = brinAdn.find("ATT")
    indexACT = brinAdn.find("ACT")
    indexATC = brinAdn.find("ATC")

    #
    while indexATT != -1 or indexACT != -1 or indexATC != -1:
        #
        if indexATT != -1:
            #
            position.append(indexATT)
            indexATT = brinAdn.find("ATT", indexATT + 3)

        if indexACT != -1:

            position.append(indexACT)
            indexACT = brinAdn.find("ACT", indexACT + 3)

        if indexATC != -1:

            position.append(indexATC)
            indexATC = brinAdn.find("ATC", indexATC + 3)
    #
    position.sort()

    #
    return position


def testTrouveFin(brinAdn):

    assert trouveFin("ACGTCGTAGCGATCGATCGTACGTCGTAGCTAG") == []
    assert trouveFin("ACTGTCGATCGTAGCTAGCGATCGATCGTAGCT") == [0]
    assert trouveFin("CGTCGATCGTAGCTAGCGATCGATCGTAGCATT") == [33]
    assert trouveFin("ACGTCGATCGTAGCGCTATTGCTAGCTAGCTAG") == [6, 17]


def trouveGene(debut: list, fin: list) -> list:
    #
    curseurDebut = 0
    curseurFin = 0
    positionGene = []

    for i in debut:
        for j in fin:
            if i < j:
                if (j-i) % 3 == 0:
                    break
        positionGene.append((i, j))

    return positionGene


def testTrouveGene(debut, fin):

    assert trouveGene([], []) == []
    assert trouveGene([0], [3]) == [(0, 3)]
    assert trouveGene([0, 6, 12], [3, 9, 15]) == [(0, 3), (6, 9), (12, 15)]
    assert trouveGene([0, 6, 10], [8, 12, 15]) == [(0, 8), (6, 12), (10, 15)]
    assert trouveGene([0, 6], [15]) == [(0, 15)]


def transcrire(brinAdn: str) -> str:

    # define a dict that maps every nucleotide to it's ARN complement
    complementArn = {
        "A": "U",
        "G": "C",
        "T": "A",
        "C": "G"
    }
    # Initialize an empty string to store the complement of the DNA strand
    brinArn = ""

    # concatenate every complementary nucleotide to generate the ARN strand
    for i in brinAdn:
        brinArn += complementArn[i]

    return brinArn


def testTranscrire(brinAdn):

    assert (transcrire("AACGTGTCTGAAGCTAGCTGGATCCTAGCGATCG") == "AAUGUGUCUGAAGCUAGCUGGAUCCUAGCGAUCG")
    assert (transcrire("TTAGCGCTAGTCTAGCTAGCTAGCTAGCTAGCTA") == "UUCGCGCUAGUCUAGCUAGCUAGCUAGCUAGCUA")
    assert (transcrire("ATCGCGATAGCGCTAGCTGATCGATCGGCTAGCT") == "AUCGCGAUAGCGCUAGCUGAUCGAUCGGCUAGCU")


def traduire(brinArn: str) -> None:

    # Initialize a cursor to track nucleotide position and an empty string that will contain the protein sequence
    curseur = 0
    protein = ""
    proteinLettre = ""
    coté = 20

    # create a string chain representing the ARN corresponding protein
    while curseur < len(brinArn) - 3:
        # from the cursor , the subsequent codon is identified
        codon = brinArn[curseur] + brinArn[curseur + 1] + brinArn[curseur + 2]
        # if it's a stop codon, no hyphen's added
        if curseur == len(brinArn) - 6:
            protein += codons_aa[codon]

        else:
            protein += codons_aa[codon] + " - "

        proteinLettre += lettreAa[codon]

        curseur += 3

    print(protein)

    indice = 0
    for i in proteinLettre:
        carre(coté, indice)
        ecrireCarre(coté, indice, i)
        indice += 1

    # saut de ligne
    sautdeligne = (indice // 15) + 2
    turtle.rt(90); turtle.fd(sautdeligne*coté)
    turtle.lt(90)


def carre(longueur: float, nombre: int) -> None:

    decalageX = (nombre % 15) * longueur
    decalageY = (nombre // 15) * longueur

    turtle.penup()

    turtle.fd(decalageX); turtle.rt(90)
    turtle.fd(decalageY); turtle.lt(90)
    turtle.pendown()

    for i in range(4):
        turtle.fd(longueur)
        turtle.rt(90)

    turtle.penup()

    turtle.rt(90); turtle.bk(decalageY)
    turtle.lt(90); turtle.bk(decalageX)


#def sauteLigne(nombre: int) -> int:
    indexDeCarre = 15 - (nombre % 15)


# def centrerGrille(longueur: int, maxLongueur: int, maxLargeur: int):
#     pass


def ecrireCarre(longueur: float, nombre: int, texte: str):

    decalageX = (nombre % 15) * longueur
    decalageY = (nombre // 15) * longueur

    turtle.fd(decalageX + longueur/2); turtle.rt(90)
    turtle.fd(decalageY + longueur/1.25); turtle.lt(90)
    turtle.pendown()

    turtle.write(texte)

    turtle.penup()
    turtle.rt(90); turtle.bk(decalageY + longueur/1.25)
    turtle.lt(90); turtle.bk(decalageX + longueur/2)


#######################################################################################################################
#                                        FONCTION PRINCIPALE ET TEST UNITAIRES
#######################################################################################################################

def adnToProtein(brinAdn):

    brinAntisens = antisens(brinAdn)

    positionGeneBrinSup = trouveGene(trouveDebut(brinAdn), trouveFin(brinAdn))
    positionGeneBrinInf = trouveGene(trouveDebut(brinAntisens), trouveFin(brinAntisens))


    genes = []

    for i in positionGeneBrinSup:
        genes.append(brinAdn[i[0]:i[1]+3])
    for i in positionGeneBrinInf:
        genes.append(brinAntisens[i[0]:i[1]+3])


    brinArn = []

    for i in genes:
        brinArn.append(transcrire(i))


    for i in brinArn:
        traduire(i)


