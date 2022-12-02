# Perl // Python
#   =~ // ==
#
# -------------
# my ($value) = @_          || указание изначальной переменной в функции
# s/A/O/g -> mAmbA => mOmbO __ /g (global find) /s (oneline value)
from icecream import ic  # библиотека для отладки
# импорт библиотек
from Bio import SeqIO
from datetime import datetime
import re
import os
import csv
import sys
import json
import configparser
import Bio.Data.CodonTable


def isProgInstalled(program):
    # vmatch2, mkvtree2, vsubseqselect2, fuzznuc, needle
    path = os.getcwd()
    for root, dirs, files in os.walk(path):
        if program in files:
            return os.path.join(root, program)


def casFinder(resDir, inFile):
    if json.loads(parametrs['useProkka'].lower()):
        repProkka = f'{resDir}\\prokka_'
        if not os.path.isdir(repProkka): os.mkdir(repProkka)
    else:
        repProkka = f'{resDir}\\prodigal_'
        if not os.path.isdir(repProkka): os.mkdir(repProkka)
    jsonCAS = ''

    nbCas = 0
    default = 0
    addToMaxSy = ''
    if parametrs['definition'].lower() == ('general' or 'g'):
        cas_db = f'{os.getcwd()}\\DEF-Class'
    elif parametrs['definition'].lower() == ('typing' or 't'):
        cas_db = f'{os.getcwd()}\\DEF-Typing'
    elif parametrs['definition'].lower() == ('subtyping' or 's'):
        cas_db = f'{os.getcwd()}\\DEF-SubTyping'

    # строка 1212 ?

    outDir_tsv = f'{resDir}\\TSV'
    if not os.path.isdir(outDir_tsv): os.mkdir(outDir_tsv)
    with open(f'{outDir_tsv}\\Cas_REPORT.tsv', 'wt') as results:
        tsv_writer = csv.writer(results, delimiter='\t')
        tsv_writer.writerow(['...', '...'])
    results.close()
    with open(f'{outDir_tsv}\\CRISPR-Cas_Systems_vicinity.tsv', 'wt') as allCas:
        tsv_writer = csv.writer(allCas, delimiter='\t')
        tsv_writer.writerow(['...', '...'])
    allCas.close()
    with open(f'{outDir_tsv}\\CRISPR-Cas_summary.tsv', 'wt') as resultsCRISPRCasSummary:
        tsv_writer = csv.writer(resultsCRISPRCasSummary, delimiter='\t')
        tsv_writer.writerow(['...', '...'])
    resultsCRISPRCasSummary.close()


def makeHTML(gff, casFile, resDir, refSeq, seqDesc, seqLen, globalAT, nbrcris, OneSpacerCris_nbr):
    pass  # заглушка ввиду открытия файла в функции


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
    # проверка параметров запуска | to_do
    function = config["Launch Function"]["Function"].split(' ')
    ic(function)
    options = {'-in': 'userfile', '-i': 'userfile',
               '-out': 'outputDirName', '-outdir': 'outputDirName',
               '-keepAll': 'keep', '-keep': 'keep',
               '-LOG': 'logOption', '-log': 'logOption',
               '-HTML': 'html', '-html': 'html',
               '-copyCSS': 'cssFile',
               '-soFile': 'so', '-so': 'so',
               '-minSeqSize': 'seqMinSize', '-mSS': 'seqMinSize',
               '-mismDRs': 'DRerrors', '-md': 'DRerrors',
               '-truncDR': 'DRtrunMism', '-t': 'DRtrunMism',
               '-minDR': 'M1', '-mr': 'M1',
               '-maxDR': 'M2', '-xr': 'M2',
               '-minSR': 'S1', '-ms': 'S1',
               '-maxSR': 'S2', '-xs': 'S2',
               '-noMism': 'mismOne', '-n': 'mismOne',
               '-percSPmin': 'Sp1', '-pm': 'Sp1',
               '-percSPmax': 'Sp2', '-px': 'Sp2',
               '-spSim': 'SpSim', '-s': 'SpSim'}  # остановился на -DBcrispr
    bool_options = ['keep', 'logOption', 'html']

    # коррекция первичных параметров согласно требованиям пользователя
    for item in function:
        if item in options:
            if options[item] not in bool_options:
                parametrs[options[item]] = function[function.index(item) + 1]
            else:
                parametrs[options[item]] = '1'
    if not isProgInstalled(parametrs['userfile']):
        sys.exit(f'Not found {parametrs["userfile"]}')
    parametrs['outputDirName'] = function[function.index("-out") + 1] if '-out' in function else 'Result'
    ic(parametrs)

    # коррекция DRs
    drTrunMism = 100 / float(parametrs['DRtrunMism'])
    drErrors = float(parametrs['DRerrors']) / 100

    # отметка времени при запуске процесса
    start_time = datetime.now().strftime("%Y-%m-%d | %H:%M:%S")
    print(f'Launch Time: {start_time}')

    # запуск работы с BIO | to_do
    # record = SeqIO.read(userfile, "fasta")
    inputfileCount = 0

    # создание папки с итогом
    outDir = f'{os.getcwd()}\\{parametrs["outputDirName"]}'
    if not os.path.isdir(outDir): os.mkdir(outDir)

    outDir_gff = f'{outDir}\\GFF'
    if not os.path.isdir(outDir_gff): os.mkdir(outDir_gff)

    # 318 - 322 ??
    # 325 - 344 для вузуализации CRISPRs and Cas genes (стоит удалить ?)
    # 348 получение размера входного файла в байтах

    # для логов
    if json.loads(parametrs['keep'].lower()):
        outDir_log = f'{outDir}\\Temporary File'
    if not os.path.isdir(outDir_log): os.mkdir(outDir_log)
    if json.loads(parametrs['logOption'].lower()) and json.loads(parametrs['keep'].lower()):
        with open(f'{outDir_log}\\log.txt', 'wt') as logfile:
            lof_write = csv.writer(logfile, delimiter='\t')
            lof_write.writerow([parametrs["userfile"], f'SIZE: {os.path.getsize(parametrs["userfile"])} bytes'])
        logfile.close()

    # JSON file 'result'
    with open(f'{outDir}\\jsonResult.json', 'wt') as jsonRes:
        json_write = csv.writer(jsonRes, delimiter='\n')
    jsonRes.close()

    casFinder(outDir, parametrs['userfile'])


# подключение базовых настроек в INI файле
config = configparser.RawConfigParser()
config.optionxform = str
config.read("settings.ini")

parametrs = {}
[parametrs.update({item[0]: item[1]}) for item in config['Base Variable'].items()]

# запуск проверки наличия программ | to_do
print(f'Welcome to {config["System Variable"]["casfinder"]}.\n')
list_programm = ['vmatch2.txt', 'mkvtree2', 'vsubseqselect2', 'fuzznuc', 'needle']
for item in list_programm:
    if isProgInstalled(item):
        print(f'_____{item} is found_____')
    # else:  # активировать, когда будут программы
        # sys.exit(f'Not found {item}')

active()  # основная ветка действий
