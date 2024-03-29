# Name: Finney Bado & Johnny Kana
# Date: 2024-03-08
# Description: This program synthesize a protein amino-acid sequence given a DNA nucleotide sequence
# The result is displayed as grid


#######################################################################################################################
#                                        DONNÉES
#######################################################################################################################

import turtle

clear(800, 600)
goto(-350, 250)

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
    "GGG": "Glycine"
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
    "GGG": "G"
}


#######################################################################################################################
#                                        SOUS-FONCTIONS ET TESTS UNITAIRES
#######################################################################################################################


def antisens(brinAdn):
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


def testAntisens():
    assert antisens("AAAAAAAAAAAAAAAAAA") == "TTTTTTTTTTTTTTTTTT"
    assert antisens("TTTTTTTTTTTTTTTTTT") == "AAAAAAAAAAAAAAAAAA"
    assert antisens("ATCGATCGATCGATC") == "GATCGATCGATCGAT"


def trouveDebut(brinAdn):
    # Initialize an empty list to store the positions of gene strand beginnings
    position = []
    # Find the position of the first occurrence of "TAC" sequence in the DNA strand (if it exists )
    index = brinAdn.find("TAC")

    # Find the positions of the next occurrences of "TAC" and add them to the list (if they exist)
    while index != -1:
        position.append(index)
        index = brinAdn.find("TAC", index + 3)

    return position


def testTrouveDebut():
    assert trouveDebut("ACGTCGTAGCGATCGATCGTAGGTCGTAGCTAG") == []
    assert trouveDebut("TACGTCGATCGTACGTAGCGATCGATCGTAGCT") == [0, 11]
    assert trouveDebut("ACGTCGATCGTATAGCTATAGCTAGCTAGGTAC") == [30]


def trouveFin(brinAdn):
    # Initialize
    position = []
    positionOrd = []
    indexATT = brinAdn.find("ATT")
    indexACT = brinAdn.find("ACT")
    indexATC = brinAdn.find("ATC")

    #
    while indexATT != -1 or indexACT != -1 or indexATC != -1:

        if indexATT != -1:
            position.append(indexATT)
            indexATT = brinAdn.find("ATT", indexATT + 3)
        else:
            pass

        if indexACT != -1:
            position.append(indexACT)
            indexACT = brinAdn.find("ACT", indexACT + 3)
        else:
            pass

        if indexATC != -1:
            position.append(indexATC)
            indexATC = brinAdn.find("ATC", indexATC + 3)
        else:
            pass

    for i in position:

        if len(positionOrd) == 0:
            positionOrd.append(i)
        else:
            index = 0
            # Find the index where the number should be inserted
            while index < len(positionOrd) and positionOrd[index] < i:
                index += 1
            # Insert the number at the correct index
            positionOrd.insert(index, i)

    return positionOrd


def testTrouveFin():
    assert trouveFin("ACGTCGTAGCGATAGACCGTACGTCGTAGCTAG") == []
    assert trouveFin("ATTGTCGATAGTAGCTAGCGATGGTTCGTAGCT") == [0]
    assert trouveFin("CGTCGATGGTAGCTAGCGAGCGAACGTAGCACT") == [30]
    assert trouveFin("ACGTCGATCGTAGCGCTATTGCTAGCTAGCTAG") == [6, 17]


def trouveGene(debut, fin):
    #
    curseurDebut = 0
    curseurFin = 0
    positionGene = []

    for i in debut:
        for j in fin:
            if i < j:
                if (j - i) % 3 == 0:
                    positionGene.append((i, j))
                    break

    return positionGene


def testTrouveGene():
    assert trouveGene([], []) == []
    assert trouveGene([0], [3]) == [(0, 3)]
    assert trouveGene([0, 6, 12], [3, 9, 15]) == [(0, 3), (6, 9), (12, 15)]
    assert trouveGene([3, 6, 10], [8, 12, 15]) == [(3, 12), (6, 12)]
    assert trouveGene([0, 6], [15]) == [(0, 15), (6, 15)]


def transcrire(brinAdn):
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


def testTranscrire():
    assert (transcrire("TACAACGTGTCTGAAGCTAGCTGGATCCTAGCGATCGACT") == "AUGUUGCACAGACUUCGAUCGACCUAGGAUCGCUAGCUGA")
    assert (transcrire("TACTTAGCGCTAGTCTAGCTAGCTAGCTAGCTAGCTAATT") == "AUGAAUCGCGAUCAGAUCGAUCGAUCGAUCGAUCGAUUAA")
    assert (transcrire("TACATGGCGATAGCGCTAGCTGATCGACCGGCTAGCTATC") == "AUGUACCGCUAUCGCGAUCGACUAGCUGGCCGAUCGAUAG")


def traduire(brinArn):
    # Initialize a cursor to track nucleotide position and an empty string that will contain the protein sequence
    curseur = 0
    protein = ""
    proteinLettre = ""
    cote = 20

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

    print(protein + "\n")

    indice = 0
    for i in proteinLettre:
        carre(cote, indice)
        ecrireCarre(cote, indice, i)
        indice += 1

    # saut de ligne
    sautdeligne = (indice // 15) + 2
    turtle.pu()
    turtle.rt(90);
    turtle.fd(sautdeligne * cote)
    turtle.lt(90)
    turtle.pd()


def carre(longueur, nombre):
    decalageX = (nombre % 15) * longueur
    decalageY = (nombre // 15) * longueur

    turtle.pu()

    turtle.fd(decalageX);
    turtle.rt(90)
    turtle.fd(decalageY);
    turtle.lt(90)
    turtle.pd()

    for i in range(4):
        turtle.fd(longueur)
        turtle.rt(90)

    turtle.pu()

    turtle.rt(90);
    turtle.bk(decalageY)
    turtle.lt(90);
    turtle.bk(decalageX)

    turtle.pd()


def ecrireCarre(longueur, nombre, texte):
    decalageX = (nombre % 15) * longueur
    decalageY = (nombre // 15) * longueur

    turtle.pu()

    turtle.fd(decalageX + longueur / 2);
    turtle.rt(90)
    turtle.fd(decalageY + longueur / 2);
    turtle.lt(90)
    turtle.pd()

    turtle.write(texte)

    turtle.pu()
    turtle.rt(90);
    turtle.bk(decalageY + longueur / 2)
    turtle.lt(90);
    turtle.bk(decalageX + longueur / 2)
    turtle.pd()


#######################################################################################################################
#                                        FONCTION PRINCIPALE ET TEST UNITAIRES
#######################################################################################################################

def adnToProtein(brinAdn):
    brinAntisens = antisens(brinAdn)

    positionGeneBrinSup = trouveGene(trouveDebut(brinAdn), trouveFin(brinAdn))
    positionGeneBrinInf = trouveGene(trouveDebut(brinAntisens), trouveFin(brinAntisens))

    genes = []

    for i in positionGeneBrinSup:
        genes.append(brinAdn[i[0]:i[1] + 3])
    for i in positionGeneBrinInf:
        genes.append(brinAntisens[i[0]:i[1] + 3])

    brinArn = []

    for i in genes:
        brinArn.append(transcrire(i))

    for i in brinArn:
        traduire(i)


#adnToProtein(adn)
