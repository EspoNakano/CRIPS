# Perl // Python
#   =~ // ==
#
# -------------
# my ($value) = @_          || указание изначальной переменной в функции
# s/A/O/g -> mAmbA => mOmbO __ /g (global find) /s (oneline value)
# -------------
# функции average, min, max, median, stdev, log2 нет смысла переписывать
from Bio import SeqIO
from tkinter import *
from tkinter.filedialog import askopenfilename
import re
import Bio.Data.CodonTable


def makeHTML(gff, casFile, resDir, refSeq, seqDesc, seqLen, globalAT, nbrcris, OneSpacerCris_nbr):
    Tk().withdraw()
    folder_selected = askopenfilename()
    print(folder_selected)
    pass


def foundInCRISPRdb(seq, start, end):
    pass  # заглушка ввиду открытия файла в функции


def repeatIDandNb(string):
    # chomp - удаление завершающей строки
    string = string.replace(' ', '')

    # далее идёт работа с открытием файла


def atpercent(string):  # calculate AT%
    if len(string) == 0:
        return 0
    else:
        res = sum([1 if n in 'AT' else 0 for _ in string])
        return res / len(string) * 100
        # return "{:.2%}".format(res / len(string))


def sequenceAlignmentMuscle(file):
    pass  # работа с файлом


def add_spacer(str_position, str_len, i, spacer):
    pass  # ???


def trans_struct_hach(position, structspacers):  # Использовать ctypes ?
    global Length, Pos1
    spacers = {}
    for i in range(structspacers):
        pos = structspacers[i]
        Pos1 = pos - position
        Length = pos
        spacers[Pos1] = pos
    return spacers


def compare_clusters(el2, el1):  # функция прмменяется только 1 раз, так что скорее всего можно встроить напрямую...
    if (el2 >= el1 - 1500) and (el2 <= el1 + 1500):
        return True
    else:
        return False


def printversion():
    pass  # надобность ?
