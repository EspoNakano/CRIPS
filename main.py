# Perl // Python
#   =~ // ==
#
# -------------
# my ($value) = @_          || указание изначальной переменной в функции
# s/A/O/g -> mAmbA => mOmbO __ /g (global find) /s (oneline value)

# импорт библиотек
from Bio import SeqIO
from tkinter import *
from datetime import datetime
from tkinter.filedialog import askopenfilename
import re
import configparser  # для чтения INI файла
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
    return 0 if len(string) == 0 else sum([1 if n in 'AT' else 0 for _ in string]) / len(string) * 100
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
    return True if (el2 >= el1 - 1500) and (el2 <= el1 + 1500) else False


def active():
    # коррекция DRs
    DRtrunMism = 100 / float(config["Base Variable"]["DRtrunMism"])
    DRerrors = float(config["Base Variable"]["DRerrors"]) / 100

    # отметка времени при запуске процесса
    start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f'Launch Time: {start_time}')

    # запуск работы с BIO | to_do
    # record = SeqIO.read(userfile, "fasta")
    inputfileCount = 0

    basename = ...  # это выходной файл ?
    outdir = ...

    ResultDir = ''


# подключение базовых настроек в INI файле
# config["Base Variable"]["SpSim"]
config = configparser.ConfigParser()
config.read("settings.ini")

print(f'Welcome to {config["System Variable"]["casfinder"]}.\n')

# проверка наличия входного файла | to_do

# запуск проверки наличия программ | to_do

# оповещение о вызове help
print('Argument -help or -h for call help')
print('For launch CasFinder use: ***\n')

# проверка параметров запуска | to_do
userfile = ...

# основная ветка действий
active()
