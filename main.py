# Perl // Python
#   =~ // ==
#
# -------------
# my ($value) = @_          || указание изначальной переменной в функции
# s/A/O/g -> mAmbA => mOmbO __ /g (global find) /s (oneline value)
from icecream import ic  # библиотека для отладки

# горячие названия файлов в PERL
# STDERR - принт ошибки в консоли
# LOG - log.txt
# JSONRES - jsonResult.json
# RESULTCCSONE - CRISPR-Cas_summary.tsv
# RESULTCCC - CRISPR-Cas_clusters.tsv

# 318 - 322 ?? -> словарь информации о повторе (шаблонизатор:?)
# 325 - 344 для вузуализации CRISPRs and Cas genes (пока игнорить)

# импорт библиотек
from Bio import SeqIO
from datetime import datetime
import re
import os
import csv
import sys
import json
import pandas as pd
import configparser
import Bio.Data.CodonTable


def isProgInstalled(program):
    # vmatch2, mkvtree2, vsubseqselect2, fuzznuc, needle
    path = os.getcwd()
    for root, dirs, files in os.walk(path):
        if program in files:
            return os.path.join(root, program)


def write_text(file, key, text):
    with open(file, key) as write_file:
        write_file.write(text)
    write_file.close()


def makeHTML(gff, casFile, resDir, refSeq, seqDesc, seqLen, globalAT, nbrcris, OneSpacerCris_nbr):
    pass  # заглушка ввиду открытия файла в функции


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

    write_text(f'{outDir_tsv}\\Cas_REPORT.tsv', 'w', '')
    write_text(f'{outDir_tsv}\\CRISPR-Cas_Systems_vicinity.tsv', 'w', '')
    write_text(f'{outDir_tsv}\\CRISPR-Cas_summary.tsv', 'a', '')

    # 1208 - 1209

    tabCCScas = []
    hashCountCas = {}

    tabCRISPRCasClusters = []
    hashCasBegin = {}
    hashCasEnd = {}
    hashCrisprBegin = {}
    hashCrisprEnd = {}

    hashCrisprDRcons = {}
    limitSide = 3000
    matchingCasType = 'CAS-TypeII'

    otherTabCas = []
    otherTabCrispr = []

    write_text(f'{outDir_tsv}\\CRISPR-Cas_clusters.tsv', 'w', '')

    # ловушка исключений

    # JSON Parse

    if json.loads(parametrs['useProkka'].lower()):
        if isProgInstalled('prokka'):
            print(f'_____{item} is found_____')
        else:  # строка 1242-1253
            if json.loads(parametrs['logOption'].lower()):
                write_text(f'{outDir_log}\\log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tprokka is not installed\n')
                write_text(f'{outDir_log}\\log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tPlease install it by following the documentation provided here: https://github.com/tseemann/prokka  OR  http://www.vicbioinformatics.com/software.prokka.shtml\n')
                write_text(f'{outDir_log}\\log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tOtherwise, please retry without Cas option (-cas)\n')
                launchCasFinder = 0
                # next ?
    elif json.loads(parametrs['useProdigal'].lower()):
        if isProgInstalled('prodigal'):
            print(f'_____{item} is found_____')
        else:  # строка 1242-1253
            if json.loads(parametrs['logOption'].lower()):
                write_text(f'{outDir_log}\\log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tprodigal is not installed\n')
                write_text(f'{outDir_log}\\log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tInstall it by following the documentation provided here: https://github.com/hyattpd/Prodigal\n')
                write_text(f'{outDir_log}\\log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tOtherwise, please retry without Cas option (-cas)\n')
                launchCasFinder = 0
                # next ?

    # STOP 1294


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


def fullReport(direct):
    pass


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
    global outDir_tsv, outDir_log

    # Входная функция
    function = config["Launch Function"]["Function"].split(' ')

    # Подготовка словаря функций
    options = {'-in': 'userfile', '-i': 'userfile',
               '-out': 'outputDirName', '-outdir': 'outputDirName',
               '-keepAll': 'keep', '-keep': 'keep',
               '-LOG': 'logOption', '-log': 'logOption',
               '-HTML': 'html', '-html': 'html',
               '-copyCSS': 'cssFile',
               '-soFile': 'so', '-so': 'so',
               '-quite': 'quiet', '-q': 'quiet',
               '-minSeqSize': 'seqMinSize', '-mSS': 'seqMinSize',
               # Detection of CRISPR arrays
               '-mismDRs': 'DRerrors', '-md': 'DRerrors',
               '-truncDR': 'DRtrunMism', '-t': 'DRtrunMism',
               '-minDR': 'M1', '-mr': 'M1',
               '-maxDR': 'M2', '-xr': 'M2',
               '-minSR': 'S1', '-ms': 'S1',
               '-maxSR': 'S2', '-xs': 'S2',
               '-noMism': 'mismOne', '-n': 'mismOne',
               '-percSPmin': 'Sp1', '-pm': 'Sp1',
               '-percSPmax': 'Sp2', '-px': 'Sp2',
               '-spSim': 'SpSim', '-s': 'SpSim',
               '-DBcrispr': 'crisprdb', '-dbc': 'crisprdb',
               '-repeats': 'repeats', '-rpts': 'repeats',
               '-DIRrepeat': 'dirRepeat', '-drpt': 'dirRepeat',
               '-flank': 'flankingRegion', '-fl': 'flankingRegion',
               '-levelMin': 'levelMin', '-lMin': 'levelMin',
               '-classifySmallArrays': 'classifySmall', '-classifySmall': 'classifySmall', '-cSA': 'classifySmall',
               '-forceDetection': 'force', '-force': 'force',
               '-fosterDRLength': 'fosteredDRLength', '-fDRL': 'fosteredDRLength',
               '-fosterDRBegin': 'fosteredDRBegin', '-fDRB': 'fosteredDRBegin',
               '-fosterDREnd': 'fosteredDREnd', '-fDRE': 'fosteredDREnd',
               '-MatchingRepeats': 'repeatsQuery', '-Mrpts': 'repeatsQuery',
               '-minNbSpacers': 'minNbSpacers', '-mNS': 'minNbSpacers',
               '-betterDetectTrunc': 'betterDetectTruncatedDR', '-bDT': 'betterDetectTruncatedDR',
               '-PercMismTrunc': 'percentageMismatchesHalfDR', '-PMT': 'percentageMismatchesHalfDR',
               # Detection of Cas clusters
               '-cas': 'launchCasFinder', '-cs': 'launchCasFinder',
               '-ccvRep': 'writeFullReport', '-ccvr': 'writeFullReport',
               '-vicinity': 'vicinity', 'vi': 'vicinity',
               '-CASFinder': 'casfinder', '-cf': 'casfinder',
               '-cpuMacSyFinder': 'cpuMacSyFinder', '-cpuM': 'cpuMacSyFinder',
               '-rcfowce': 'rcfowce',
               '-definition': 'definition', '-def': 'definition',
               '-gffAnnot': 'userGFF', '-gff': 'userGFF',
               '-proteome': 'userFAA', '-faa': 'userFAA',
               '-cluster': '', '-ccc': 'clusteringThreshold',
               '-getSummaryCasfinder': 'gscf', '-gscf': 'gscf',
               '-geneticCode': 'genCode', '-gcode': 'genCode',
               '-metagenome': 'metagenome', '-meta': 'metagenome',
               # [Use Prokka instead of Prodigal (default option)]
               # Prokka (https://github.com/tseemann/prokka) must be installed in order to use following options
               '-useProkka': 'useProkka', '-prokka': 'useProkka',
               '-cpuProkka': 'cpuProkka', '-cpuP': 'cpuProkka',
               '-ArchaCas': 'kingdom', '-ac': 'kingdom'}
    bool_options = ['keep', 'logOption', 'html', 'kingdom']

    # Коррекция первичных параметров согласно требованиям пользователя
    for item in function:
        if item in options:
            if options[item] not in bool_options:
                parametrs[options[item]] = function[function.index(item) + 1]
            else:
                if item == ('-ArchaCas' or '-ac'):
                    parametrs[options[item]] = 'Archaea'
                else:
                    parametrs[options[item]] = '1'
        else:
            if item[0] == '-':
                sys.exit(f'Error argument {item}')
    parametrs['outputDirName'] = function[function.index("-out") + 1] if '-out' in function else 'Result'
    if json.loads(parametrs['force'].lower()):
        parametrs['M1'] = parametrs['fosteredDRLength']
    if json.loads(parametrs['useProkka'].lower()):
        parametrs['useProdigal'] = '0'
    else:
        parametrs['useProdigal'] = '1'
    if parametrs['kingdom'] == 'Archaea':
        parametrs['launchCasFinder'] = '1'

    # Проверка наличия входного файла
    if not isProgInstalled(parametrs['userfile']):
        sys.exit(f'Not found {parametrs["userfile"]}')

    # коррекция DRs
    drTrunMism = 100 / float(parametrs['DRtrunMism'])
    drErrors = float(parametrs['DRerrors']) / 100

    # отметка времени при запуске процесса
    start_time = datetime.now().strftime("%d-%m-%Y | %H:%M:%S")
    print(f'Launch Time: {start_time}')

    # запуск работы с BIO | to_do
    # record = SeqIO.read(userfile, "fasta")
    inputfileCount = 0

    # создание папки с итогом
    outDir = f'{os.getcwd()}\\{parametrs["outputDirName"]}'
    outDir_tsv = f'{outDir}\\TSV'
    if not os.path.isdir(outDir): os.mkdir(outDir)
    if not os.path.isdir(outDir_tsv): os.mkdir(outDir_tsv)

    outDir_gff = f'{outDir}\\GFF'
    if not os.path.isdir(outDir_gff): os.mkdir(outDir_gff)

    # для логов (ПОСМОТРЕТЬ LOGGING)
    if json.loads(parametrs['keep'].lower()):
        outDir_log = f'{outDir}\\Temporary File'
    if not os.path.isdir(outDir_log): os.mkdir(outDir_log)
    if json.loads(parametrs['logOption'].lower()) and json.loads(parametrs['keep'].lower()):
        with open(f'{outDir_log}\\log.txt', 'w') as logfile:
            lof_write = csv.writer(logfile, delimiter='\t')
            lof_write.writerow([parametrs["userfile"], f'SIZE: {os.path.getsize(parametrs["userfile"])} bytes'])
        logfile.close()

    # JSON file 'result'
    write_text(f'{outDir}\\jsonResult.json', 'w', '{\n')
    write_text(f'{outDir}\\jsonResult.json', 'a', f'Date:  {start_time}  [dd-mm-yy | hh-mm-ss]\nVersion:  python\n')
    write_text(f'{outDir}\\jsonResult.json', 'a', f'Command:  {config["Launch Function"]["Function"]}\n')
    write_text(f'{outDir}\\jsonResult.json', 'a', 'Sequences:\n[{')

    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}\\log.txt', 'a',
                   f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tResults will be stored in {outDir}\n')

    if not json.loads(parametrs['quiet'].lower()):
        print(f'{outDir_log}\\log.txt', 'a',
              f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tResults will be stored in {outDir}\n')


    allFoundCrisprs = 0
    allCrisprs = 0
    nbrAllCas = 0

    actualMetaOptionValue = parametrs['metagenome']

    # разобраться с этими принтами
    write_text(f'{outDir_tsv}\\CRISPR-Cas_summary.tsv', 'w',
               'Sequence(s)\tCRISPR array(s)\tNb CRISPRs\tEvidence-levels\tCas cluster(s)\tNb Cas\tCas Types/Subtypes\n')
    if json.loads(parametrs['clusteringThreshold'].lower()) and json.loads(parametrs['launchCasFinder'].lower()):
        write_text(f'{outDir_tsv}\\CRISPR-Cas_clusters.tsv', 'w',
                   'Sequence\tCluster Id\tNb CRISPRs\tNb Cas\tStart\tEnd\tDescription\n')

    # подвязать pandas к созданию XLS
    # if json.loads(parametrs['classifySmall'].lower()):
    #    with open(f'{outDir}\\smallArraysReclassification.xls', 'wt') as smallArrays:
    #        smallArrays.write('text')
    #    smallArrays.close()

    # цикл 433-712 | работает с BIO
    casFinder(outDir, parametrs['userfile'])
    # 679 - 702
    analyzedSequences = f'{outDir}\\analyzedSequences'
    write_text(f'{outDir}\\jsonResult.json', 'a', ']\n')
    write_text(f'{outDir}\\jsonResult.json', 'a', '}\n')

    if json.loads(parametrs['launchCasFinder'].lower()) and json.loads(parametrs['writeFullReport'].lower()):
        fullReport(outDir)

    # 735 - 885


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
