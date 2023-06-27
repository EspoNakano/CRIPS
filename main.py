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
import re
import os
import csv
import sys
import json
import math
import shlex
import subprocess
import pandas as pd
import configparser
import Bio.Data.CodonTable

from Bio import SeqIO
from pathlib import Path
from io import BufferedReader
from datetime import datetime


def isProgInstalled(program):
    PATH = os.getenv('PATH').split(':')  # Get user's Path
    found = False
    for path in PATH:
        if os.access(os.path.join(path, program), os.X_OK):
            found = True
            break
    if program.startswith('muscle'):
        if not found:
            print(f"\nThe program {program} cannot be found on your system\n")
            print(f"Have you installed it? Your PATH variable contains: {os.getenv('PATH')}\n\n")
            print(f"\nPlease install {program}\n")
            exit(EX_CONFIG)
    return found


def atpercent(string):  # calculate AT%
    return 0 if len(string) == 0 else sum([1 if n in 'AT' else 0 for n in string]) / len(string) * 100


def callmkvtree(inputfile, indexname):
    if not checkdbfile(inputfile, f"{indexname}.prj"):
        system_command = f"mkvtree2 -db {inputfile} -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1"
        makesystemcall(system_command)
        if json.loads(parametrs['logOption'].lower()):
            write_text(f'{outDir_log}/log.txt', 'a',
                       f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\t{system_command}\n')


def checkdbfile(inputfile, prjfile):
    dbfile = ''
    if os.path.exists(prjfile):
        try:
            with open(prjfile, "r") as PRJFILEPTR:
                for line in PRJFILEPTR:
                    match = re.search(r'^dbfile=(\S+) (\d+)', line)
                    if match:
                        if dbfile == '':
                            dbfile = match.group(1)
                            dbfilesize = int(match.group(2))
                            if dbfile == inputfile:
                                try:
                                    size = os.path.getsize(dbfile)
                                    if size == dbfilesize:
                                        return 1
                                except OSError:
                                    pass
            return 0
        except IOError:
            return 0
    else:
        return 0


def makesystemcall(arg_string) -> None:
    args = arg_string.split(' ')
    retcode = subprocess.run(args=args)
    if retcode.returncode != 0:
        print(f"failure: \"{arg_string}\", error code {retcode}", file=sys.stderr)
        sys.exit(1)


def trans_data(file):
    elem_nbr = 0
    repetitions = []
    with open(file) as fd:
        for l in fd:
            if l.strip() == "":
                continue
            line = l.strip()
            line = line.replace(">", "> ")
            temp = line.split()
            if temp[0].startswith(">"):
                rep = Rep()
                rep.Length = temp[1]
                rep.Pos1 = temp[2]
                repetitions.append(rep)
            else:
                repetitions[elem_nbr].DRseq = temp[0]
                elem_nbr += 1
    fd.close()
    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}/log.txt', 'a',
                   f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\t'
                   f'Getting results from vmatch and transform vmatch output file...\n')
    return repetitions


def write_text(file, key, text):
    with open(file, key) as write_file:
        write_file.write(text)
    write_file.close()


def extractsequence(seqencesdefinition):  # *seqencesdefinition
    for count in range(0, len(seqencesdefinition), 2):
        if seqencesdefinition[count] - 500 > 0:
            seqencesdefinition[count] -= 500
        else:
            seqencesdefinition[count] = 0
        if seqencesdefinition[count + 1] + 500 >= seq.length():
            seqencesdefinition[count + 1] = seq.length() - 1
        else:
            seqencesdefinition[count + 1] += 500
        number = (count + 2) / 2
        seqfile = f'seq_v{number}'
        subselectoptions = ["-range", seqencesdefinition[count], seqencesdefinition[count + 1], inputfile, ">", seqfile]
        makesystemcall("vsubseqselect2 " + " ".join(subselectoptions))
    res = [*seqencesdefinition]
    return res


def casFinder(ResultDir, inputfile, seqDesc, RefSeq, nbrcris, kingdom, casfinder_opt):
    if json.loads(parametrs['useProkka'].lower()):
        repProkka = f'{resDir}/prokka_'
        if not os.path.isdir(repProkka): os.mkdir(repProkka)
    else:
        repProkka = f'{resDir}/prodigal_'
        if not os.path.isdir(repProkka): os.mkdir(repProkka)
    json = f'{resDir}/casfinder_/{RefSeq}'
    jsonCAS = ''

    nbCas = 0
    default = 0
    addToMacSy = ''
    if parametrs['definition'].lower() == ('general' or 'g'):
        cas_db = f'{os.getcwd()}/DEF-Class'
        addToMacSy += '--min-genes-required General-Class1 1 '
        addToMacSy += '--min-genes-required General-Class2 1 '
    elif parametrs['definition'].lower() == ('typing' or 't'):
        cas_db = f'{os.getcwd()}/DEF-Typing'
        # addToMacSy += '--min-genes-required CAS 1 --min-genes-required CAS-TypeI 1 --min-genes-required CAS-TypeII 1
        # --min-genes-required CAS-TypeIII 1 --min-genes-required CAS-TypeIV 1 --min-genes-required CAS-TypeV 1 '
    elif parametrs['definition'].lower() == ('subtyping' or 's'):
        cas_db = f'{os.getcwd()}/DEF-SubTyping'

    results = f'{outDir_tsv}/Cas_REPORT.tsv'
    allCas = f'{outDir_tsv}/CRISPR-Cas_Systems_vicinity.tsv'

    write_text(results, 'w', '')
    write_text(allCas, 'w', '')
    write_text(f'{outDir_tsv}/CRISPR-Cas_summary.tsv', 'a', '')

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

    write_text(f'{outDir_tsv}/CRISPR-Cas_clusters.tsv', 'w', '')

    # ловушка исключений

    # JSON Parse библиотеки

    if json.loads(parametrs['useProkka'].lower()):
        if isProgInstalled('prokka'):
            if json.loads(parametrs['quiet'].lower()):
                print(f'_____{item} is found_____')
        else:  # строка 1242-1253
            if json.loads(parametrs['logOption'].lower()):
                write_text(f'{outDir_log}/log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tprokka is not installed\n')
                write_text(f'{outDir_log}/log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tPlease install it by following the '
                           f'documentation provided here: https://github.com/tseemann/prokka  OR  '
                           f'http://www.vicbioinformatics.com/software.prokka.shtml\n')
                write_text(f'{outDir_log}/log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tOtherwise, please retry without Cas '
                           f'option (-cas)\n')
                launchCasFinder = 0
                # next ?
    elif json.loads(parametrs['useProdigal'].lower()):
        if isProgInstalled('prodigal'):
            if json.loads(parametrs['quiet'].lower()):
                print(f'_____{item} is found_____')
        else:  # строка 1242-1253
            if json.loads(parametrs['logOption'].lower()):
                write_text(f'{outDir_log}/log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tprodigal is not installed\n')
                write_text(f'{outDir_log}/log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tInstall it by following the '
                           f'documentation provided here: https://github.com/hyattpd/Prodigal\n')
                write_text(f'{outDir_log}/log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tOtherwise, please retry without '
                           f'Cas option (-cas)\n')
                launchCasFinder = 0
                # next ?
    if isProgInstalled('macsyfinder'):
        if json.loads(parametrs['quiet'].lower()):
            print('macsyfinder installation is...........OK \n')
        else:
            print('\n_\nmacsyfinder is not installed .........\n___\n')
            print('Please install it by following the documentation provided here: '
                  'https://github.com/gem-pasteur/macsyfinder\n')
            print('Otherwise, please retry without Cas option (-cas)\n\n')
            launchCasFinder = 0
            # next ?

    options = ''
    if json.loads(parametrs['quiet'].lower()):
        options += '--quiet '
    if json.loads(parametrs['fast'].lower()):
        options += '--fast --rawproduct --norrna --notrna '
    if json.loads(parametrs['metagenome'].lower()):
        options += '--metagenome '
    prokka = f'prokka {options}'
    prokka += f'--cpus {parametrs["cpuProkka"]} --kingdom {parametrs["kingdom"]} --gcode {parametrs["genCode"]} ' \
              f'--outdir {parametrs["outputDirName"]} --prefix {parametrs["RefSeq"]} {parametrs["userfile"]}'

    # 1345-1390

    json = f'{casDir}/results.macsyfinder.json'

    # 1392-1400

    resultsCRISPRs = f'{outDir_log}/Crisprs_REPORT.tsv'
    resultsTmpCRISPRs = f'{outDir_log}/TmpCrisprs_REPORT.tsv'

    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}/log.txt', 'a',
                   f'{datetime.now().strftime("%H:%M:%S")}\tWriting in file Cas_REPORT.tsv\n')

    write_text(results, 'a', '############################################\n')
    write_text(results, 'a', f'{RefSeq} ({seqDesc})\n')
    write_text(results, 'a', '--------------------------------------------\n')


def find_clusters(*repetitions):
    tabseq = []
    count = 0
    nb_clust = 0  # nb_clust = 2*nbr_clusters

    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}/log.txt', 'a',
                   f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\t'
                   f'Find clusters and store their start and end positions\n')

    tabseq.append(repetitions[count].Pos1)  # start position of the first cluster
    while count < len(repetitions) - 1:
        if not compare_clusters(repetitions[count].Pos1, repetitions[count + 1].Pos1):
            # a new cluster is found
            if repetitions[count].Pos1 != tabseq[nb_clust]:
                nb_clust += 1
                tabseq.append(repetitions[count].Pos1)
            else:
                nb_clust += 1
                tabseq.append(repetitions[count].Pos1 + 60)
            nb_clust += 1
            tabseq.append(repetitions[count + 1].Pos1)
        count += 1
    if len(tabseq) % 2 == 0:
        if repetitions[count].Pos1 != tabseq[nb_clust]:
            nb_clust += 1
            tabseq.append(repetitions[count].Pos1)
        else:
            nb_clust += 1
            tabseq.append(repetitions[count].Pos1 + 60)
    return tabseq


def DR_occ_rev(DR, rep):
    from Bio.Seq import Seq
    DRocc = 0
    seqDR = Seq(DR)
    revDR = seqDR.reverse_complement()
    for i in rep:
        if (DR == i.DRseq) or (revDR == i.DRseq):
            DRocc += 1
    return DRocc


def check_DR(DR, DR_cand):
    i = 0
    stop = 0
    while i <= len(DR_cand) - 1 and stop == 0:
        if DR == DR_cand[i]:
            stop = 1
        i += 1
    return stop


# def write_clusters(RefSeq, rep):
#     tabseq = find_clusters(rep)
#
#     if json.loads(parametrs['logOption'].lower()):
#         write_text(f'{outDir_log}/log.txt', 'a',
#                    f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tFind CRISPRs candidates and check DRs...\n')
#
#     count, seqbeg, seqend, i, j, DR, DRocc, nbrcris, DRlength = 0, 0, 0, 0, 0, None, None, 0, 0
#     i = 0
#     j = 0
#     nbrcris = 0
#     OneSpacerCris_nbr = 0
#     modf = -1
#     modifcode = []
#
#     tabseq = extractsequence(tabseq)


def compute_Crispr_score(crisprfile, DRlength):
    score = 0
    lines = []
    temp = []
    penal = 0
    TotERR = 0
    nb = 0
    troncated = 0

    with open(crisprfile, 'r') as fd:
        lines = fd.readlines()

    count = 19
    if len(lines) <= 5:
        score = 100000

    for count in range(12, len(lines) - 1):
        nb += 1
        if lines[count].startswith('#'):
            continue
        if lines[count].isspace():
            continue
        lines[count] = lines[count].strip()
        if lines[count].startswith('Start'):
            continue
        temp = lines[count].split()

        if json.loads(parametrs['force'].lower()):
            if temp[-1].startswith('G'):
                score -= 2
            if len(temp[-1]) >= 30:
                score -= 2
            if temp[-1].endswith('AA'):
                score -= 2

        if temp[-2] != '.':
            TotERR += int(temp[-2])
            penal = 1 + int(temp[-2]) / DRlength
            score += penal
            if count == 19:
                troncated = penal
            if count == len(lines) - 4:
                if penal > troncated:
                    troncated = penal
    score -= troncated
    TotERR /= (nb * DRlength)
    return score, TotERR


def Find_theCrispr(index_name, DR, count, seq_beg, seq_end, crispr_file):
    DR_length = len(DR)
    sim_DRs, ref_fals_spacers, spacers = definespacers(DR_length, crispr_file, index_name)
    Fals_spacers = ref_fals_spacers
    if len(Fals_spacers) >= 0:
        for h in range(len(Fals_spacers)):
            Fals_spacers[h] += seq_beg
    spacers_H = trans_struct_hach(seq_beg, spacers)
    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}/log.txt', 'a',
                   f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tCheck short CRISPRs arrays and check spacers '
                   f'alignment using Muscle....\n')
    good_crispr = 0
    o = len(spacers) - 1
    if min_nb_spacers <= 0:
        min_nb_spacers = 1
    if len(spacers) >= (min_nb_spacers - 1):
        if len(spacers) == 0:
            good_crispr = check_short_crispr(DR, f"spacers{index_name}")
        else:
            good_crispr = checkspacersAlignMuscle(f"spacers{index_name}")
        if good_crispr:
            Cris_beg = seq_beg + spacers[0].pos1 - DR_length
            Cris_end = seq_beg + spacers[-1].pos1 + spacers[-1].length + DR_length - 1
    r = len(Fals_spacers) - 1
    return good_crispr, sim_DRs, Cris_beg, Cris_end, spacers_H, len(spacers), Fals_spacers


def write_clusters(RefSeq, *rep):
    tabseq = find_clusters(*rep)
    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}/log.txt', 'a',
                   f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tFind CRISPRs candidates and check DRs...\n')
    i = 0
    j = 0
    nbrcris = 0
    OneSpacerCris_nbr = 0
    modf = -1
    modifcode = []
    tabseq = extractsequence(tabseq)
    count = 1
    while count <= len(tabseq) // 2 + 1:
        DR_cand = []
        DR_cand_occ = []
        seqbeg = tabseq[2 * count - 2]
        seqend = tabseq[2 * count - 1]
        DR = rep[i].DRseq
        DR_cand.append(DR)
        occ_DR_max = DR_occ_rev(DR, *rep)
        DR_cand_occ.append(DR)
        i += 1
        while i <= len(rep) and rep[i].Pos1 <= seqend:
            DR = rep[i].DRseq
            if not check_DR(DR, *DR_cand):
                DR_cand.append(DR)
                DRocc = DR_occ_rev(DR, *rep)
                if occ_DR_max < DRocc:
                    occ_DR_max = DRocc
                    DR_cand_occ = []
                    DR_cand_occ.append(DR)
                else:
                    if occ_DR_max == DRocc:
                        if DR != DR_cand_occ[0]:
                            DR_cand_occ.append(DR)
            i += 1
        count += 1

        indexname = "seq_v" + str(count)
        crisprfile = "crispr_result_" + str(count)
        actual_path = os.getcwd()
        bestscore = 100000
        bestindex = 0
        c = 0
        TotDRerr = 0
        TotDRerr_best = 0

        while c <= len(DR_cand_occ) - 1:
            score = None

            if len(DR_cand_occ) > 1:
                crisprfile = "crispr_result_" + str(count) + "_" + str(c)

            DRlength = len(DR_cand_occ[c])
            err = 0
            if betterDetectTruncatedDR:
                charSeqDR = list(DR_cand_occ[c])
                newDRpattern = ''
                for k in range(DRlength):
                    if k < int(DRlength / 2):
                        newDRpattern += charSeqDR[k]
                    else:
                        newDRpattern += "n"
                err = betterDetectTruncatedDR
                fuzznucoptions = ["-sequence", indexname, "-pattern", newDRpattern, "-pmismatch", err, "-outfile",
                                  crisprfile]
            else:
                err = int(DRlength / DRtrunMism)
                fuzznucoptions = ["-sequence", indexname, "-pattern", DR_cand_occ[c], "-pmismatch", err, "-outfile",
                                  crisprfile]
            fuzznucoptions.append("-auto")
            makesystemcall("fuzznuc " + " ".join(fuzznucoptions))

            if json.loads(parametrs['logOption'].lower()):
                write_text(f'{outDir_log}/log.txt', 'a',
                           f'\n{datetime.now().strftime("%H:%M:%S")}\t mkvtree2 -db {inputfile} '
                           f'fuzznuc {fuzznucoptions} \n')

            if len(DR_cand_occ) > 0:
                score, TotDRerr = compute_Crispr_score(crisprfile, len(DR_cand_occ[c]))
                if score <= bestscore:
                    bestscore = score
                    bestindex = c
                    TotDRerr_best = TotDRerr
            else:
                bestscore, TotDRerr_best = compute_Crispr_score(crisprfile, len(DR_cand_occ[0]))
            c += 1

        if len(DR_cand_occ) > 0:
            crisprfile = f'crispr_result_{count}'
            if len(DR_cand_occ) > 0:
                rename(f'{crisprfile}_{bestindex}', crisprfile)
            DR = DR_cand_occ[bestindex]
        else:
            DR = DR_cand_occ[0]

        (criskOK, simDRs, CrisprBeg, CrisprEnd, RefSpacersH, nbspacers,
         RefFalsSpacers) = Find_theCrispr(indexname, DR, count, seqbeg, seqend, crisprfile)

        if criskOK:
            if TotDRerr_best > DRerrors:
                criskOK = False  # DC - replace 0.2 by DRerrors (assuming DRerrors is a variable)
        nbrcris = nbrcris + criskOK

        if criskOK:
            Crispr_file = fill_in_crisprfile(simDRs, RefSeq, ResultDir, nbrcris, CrisprBeg, CrisprEnd, DR,
                                             nbspacers, "spacers" + indexname, RefFalsSpacers, RefSpacersH)

            FalsSpacers = RefFalsSpacers
            if len(FalsSpacers) != -1:
                modf = len(FalsSpacers)
            if nbspacers <= 1:
                OneSpacerCris_nbr += 1  # DC replaced '$nbspacers <= 1 || $simDRs == 0' by '$nbspacers <= 3' or '$nbspacers <= 2'

            actual_path_after_fill_in = os.getcwd()

            if json.loads(parametrs['logOption'].lower()):
                write_text(f'{outDir_log}/log.txt', 'a',
                           f'\n{datetime.now().strftime("%H:%M:%S")}\t Nb spacers = {nbspacers} , '
                           f'similarity DRs = {simDRs} , Hypotheticals = {OneSpacerCris_nbr}\n')
            if modf != -1:
                if json.loads(parametrs['logOption'].lower()):
                    write_text(f'{outDir_log}/log.txt', 'a', f'\n{datetime.now().strftime("%H:%M:%S")}\t '
                                                             f'Actual path directory before CRISPR modification: '
                                                             f'{actual_path} ...\n')
                    write_text(f'{outDir_log}/log.txt', 'a', f'\n{datetime.now().strftime("%H:%M:%S")}\t '
                                                             f'Actual path directory before CRISPR modification: '
                                                             f'{actual_path} ...\n')

                os.chdir(RefSeq)

                # DC - 07/2017 - End of DANGEROUS MODIFICATION

                hyp_cris_nbr = OneSpacerCris_nbr
                crisprs_nbr = nbrcris
                # --------------------------------------

                while modf >= 0:

                    if json.loads(parametrs['logOption'].lower()):
                        write_text(f'{outDir_log}/log.txt', 'a', f'\n{datetime.now().strftime("%H:%M:%S")}\t '
                                                                 f'Modification would be performed (Nb good CRISPRs = '
                                                                 f'{crisprs_nbr}; Nb hypothetical = {hyp_cris_nbr})'
                                                                 f'...\n)')
                    if crisprs_nbr - hyp_cris_nbr >= 1:
                        s = 1
                        while s <= crisprs_nbr:
                            # DC - 05/2017 - replace $inputfile by $inputfileTmp (in the whole 'while' loop)
                            inputfileTmp = f"{RefSeq}_Crispr_{s}"  # DC
                            if os.path.exists(inputfileTmp):
                                s, crisprs_nbr, OneSpacerCris_nbr = modify_files(inputfileTmp, ResultDir,
                                                                                 RefSeq, s, crisprs_nbr,
                                                                                 OneSpacerCris_nbr, modf)
                                nbrcris = crisprs_nbr
                            s += 1
                    modf -= 1
            modf -= 1
            if json.loads(parametrs['logOption'].lower()):
                write_text(f'{outDir_log}/log.txt', 'a', f'\n{datetime.now().strftime("%H:%M:%S")}\t Actual path '
                                                         f'directory after CRISPR modification: $actual_path ...\n')
            os.chdir(actual_path_after_fill_in)  # DC - 07/2017 - IMPORTANT addition DC removed chdir($actual_path);
    return nbrcris, OneSpacerCris_nbr


def foundInCRISPRdb(seq, start, end):
    pass  # заглушка ввиду открытия файла в функции


def modify_files(inputfile, ResultDir, RefSeq, rank, crisprs_nbr, HYPCris_nbr, modf):
    actual_pathDC = os.getcwd() # DC
    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}/log.txt', 'a', f'\n{datetime.now().strftime("%H:%M:%S")}\t Modify CRISPRs files '
                                                 f'generated ({inputfile})...........\nActual path: $actual_pathDC\n')
    dir = os.path.join(ResultDir, RefSeq)
    os.chdir(dir)

    with open(inputfile) as f:
        lines = f.readlines()

    newCrisNbr = crisprs_nbr
    temp = lines[15].split(":")
    temp2 = temp[1].split()
    begPos = temp2[1]
    endPos = temp[2]
    temp = lines[16].split(":")
    temp2 = temp[2].split()
    DRlength = temp2[1]
    nb_spacers = temp[3]
    j = 20 + 2 * int(nb_spacers)

    if '#' in lines[j]:
        pass
    else:
        temp = lines[j].split(":")
        nb_div = temp[1].split()
        if modf != -1:
            divi = temp[2:]
            divi.pop()
            newCrisNbr += len(nb_div)
        line = find_Sp_lines(j,inputfile,nb_div[0])
        (file_new, file_old, l1, l2, id, CrisprBeg, CrisprEnd, nbspacers, Spfile_old, Spfile_new) = [None] * 10
        file_old = inputfile
        Spfile_old = f"Spacers_{rank}"
    (nbsp1, nbsp2) = (None, None)

    l1 = 20
    l2 = line
    id = rank
    Spfile_new = "Spacers_test_" + str(id)
    CrisprBeg = begPos
    CrisprEnd = nb_div[0]
    CrisprEnd = CrisprEnd + "\n"
    nbsp1 = (l2 - l1) / 2

    if nbsp1 >= 1:
        file_new = inputfile + "test" + str(id)
        crisprs_nbr = create_file(file_new, file_old, l1, l2, id, CrisprBeg, CrisprEnd, nbsp1, modf, 1, crisprs_nbr, divi)
        create_spFile(dir, ResultDir, RefSeq, Spfile_new, Spfile_old, sp_prev, nbsp2)
    else:
        crisprs_nbr -= 1

    rank = id
    if nbsp1 >= 1 and nbsp2 >= 1:
        for k in range(crisprs_nbr - 1, rank - 1, -1):
            l = k + 1
            if os.path.exists(f"{RefSeq}Crispr{k}"):
                os.rename(f"{RefSeq}Crispr{k}", f"{RefSeq}Crispr{l}")
                os.rename(f"Spacers_{k}", f"Spacers_{l}")
    else:
        for k in range(crisprs_nbr - 1, rank - 1, -1):
            l = k - 1
            if os.path.exists(f"{RefSeq}Crispr{k}"):
                os.rename(f"{RefSeq}Crispr{k}", f"{RefSeq}Crispr{l}")
                os.rename(f"Spacers_{k}", f"Spacers_{l}")
    l = rank - 1

    if nbsp2 >= 1:
        if nbsp1 >= 1:
            l = rank
        os.unlink(f"{RefSeq}_Crispr_{rank}")
        if os.path.exists(f"{inputfile}_test_{rank}"):
            os.rename(f"{inputfile}_test_{rank}", f"{RefSeq}_Crispr_{l}")
        # ___spacers files ----------
        Spfile_new = f"Spacers_test_{rank}"
        Spfile = f"Spacers_{l}"
        os.rename(Spfile_new, Spfile)
    l = rank - 1

    if nbsp1 >= 1:
        os.unlink(f"{RefSeq}_Crispr_{l}")
        if os.path.exists(f"{inputfile}_test_{l}"):
            os.rename(f"{inputfile}_test_{l}", f"{RefSeq}_Crispr_{l}")
        # ___spacers files ----------
        Spfile_new = f"Spacers_test_{l}"
        Spfile = f"Spacers_{l}"
        os.rename(Spfile_new, Spfile)

    return rank, crisprs_nbr, HYPCris_nbr


def repeatIDandNb(repeat):
    repeat = repeat.strip().replace(" ", "")
    file = "Repeat_List.csv"
    line = ""
    id = ""
    number = 0
    orientation = ""
    hashID = {}
    hashNB = {}
    hashOrientation = {}
    if os.path.exists(file):
        with open(file, 'r') as f:
            for line in f:
                line = line.strip().replace(" ", "")
                temp, id, number, orientation = line.split(";")
                hashID[temp] = id
                hashNB[temp] = number
                hashOrientation[temp] = orientation
    if repeat in hashID:
        id = hashID[repeat]
        number = int(hashNB[repeat])
        orientation = hashOrientation[repeat]
    else:
        id = "Unknown"
        number = 0
        orientation = ""
    return id, number, orientation


def sequenceAlignmentMuscle(file):
    pass  # работа с файлом


def add_spacer(pos, length, i, *spacers):
    pos = pos + 1
    spacers[i] = Rep()
    spacers[i].Pos1(pos)
    spacers[i].Length(length)
    return spacers


def fullReport(direct):
    pass


def entropy(file):
    with open(file) as f:
        lines = f.readlines()

    words = {}
    total = 0
    text = []
    seqLength = 0
    tableRows = []
    tableCols = []

    table = []
    countLine = 0

    for line in lines:
        line = line.strip()
        words = re.split('[^a-zA-Z]+', line)
        if '>' not in line:
            seqLength = len(line)
            for i in range(seqLength):
                value = line[i]
                table[countLine][i] = value

                # Instantiate tableCols with values from columns
                tableCols[i][countLine] = value

            countLine += 1
    sumEntropy = 0
    for j in range(seqLength):
        sizeTable = len(tableCols[j])
        element = {}
        for elem in tableCols[j]:
            if elem == "-":
                element[elem] = 0
            else:
                if elem in element:
                    element[elem] += 1
                else:
                    element[elem] = 1
        entropy = 0
        for word in tableCols[j]:
            if word != "-":
                prob = element[word] / sizeTable
                entropy += math.log2(prob)
        entropy *= -1
        entropy = entropy / sizeTable
        sumEntropy = sumEntropy + (1 - entropy)
    finalResult = (sumEntropy / seqLength) * 100

    if finalResult < 0:
        return 0
    else:
        return finalResult


def fastaAlignmentMuscleOne(file_path):
    result = "fastaMuscle_" + file_path
    try:
        prog = isProgInstalled('muscle')
        muscle = f"muscle -in {file_path} -out {result}"
        if json.loads(parametrs['quiet'].lower()):
            muscle += " -quiet"
        if prog:
            makesystemcall(muscle)
    except:
        print("An error occurred in function fastaAlignmentMuscle")
    return result


def repeatDirection(id):
    id = id.strip().replace(" ", "")
    file = "repeatDirection.txt"
    line = ""
    orientation = ""
    hashID = {}
    if os.path.exists(file):
        with open(file, 'r') as f:
            for line in f:
                line = line.strip()
                temp, orientation = line.split("\t")
                temp = temp.strip().replace(" ", "")
                hashID[temp] = orientation
    if id in hashID:
        return hashID[id]
    else:
        return "Unknown"


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


def compute_Crisprscore(crisprfile, DRlength):
    score = 0
    lines = []
    temp = []
    penal = 0
    TotERR = 0
    nb = 0
    truncated = 0
    TotERR = 0
    nb = 0
    with open(crisprfile) as fd:
        lines = fd.readlines()
    count = 19
    if len(lines) <= 5:
        score = 100000
    count = 12
    while count < len(lines) - 1:
        nb += 1
        if re.match(r'^#', lines[count]):
            count += 1
            continue
        if re.match(r'^\s+$', lines[count]):
            count += 1
            continue
        lines[count] = lines[count].strip()
        if re.match(r'^Start', lines[count]):
            count += 1
            continue
        temp = re.split(r'\s+', lines[count])

        if force:
            if re.match(r'^G', temp[-1]):
                score -= 2
            if len(temp[-1]) >= 30:
                score -= 2
            if re.match(r'AA.$', temp[-1]):
                score -= 2

        if temp[-2] != ".":
            TotERR += temp[-2]
            penal = 1 + temp[-2] / DRlength
            score += penal
            if count == 19:
                troncated = penal
            if count == len(lines) - 4:
                if penal > troncated:
                    troncated = penal
        count += 1
    score -= troncated
    TotERR /= (nb * DRlength)
    return score, TotERR


def fastaAlignmentMuscle(file):
    result = file + "_fasta"
    try:
        prog = isProgInstalled("muscle")
        muscle = "muscle -in " + file + " -out " + result
        if json.loads(parametrs['quiet'].lower()):
            muscle += " -quiet"
        if prog:
            makesystemcall(muscle)
            if json.loads(parametrs['logOption'].lower()):
                write_text(f'{outDir_log}/log.txt', 'a', f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\t'
                                                         f'{muscle}\n')
    except:
        print("An error occurred in function fastaAlignmentMuscle")
    return result


def getSp(file_old, l1):
    temp = []
    with open(file_old, 'r') as fd:
        for line_num, line in enumerate(fd):
            if line_num == l1:
                temp = line.strip().split()
    sp = temp[1]
    return sp


def copy_spacers(SpacersFile, dir_path, id):
    dir_path = os.path.join(dir_path, f"Spacers_{id}")
    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}/log.txt', 'a', f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\t'
                                                 f'Copy Spacers files in {dir}...........\n')
    cop = f"cp -R {SpacersFile} {dir_path}"
    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}/log.txt', 'a', f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\t'
                                                 f'{cop}...........\n')
    os.system(cop)


def find_Sp_lines(j, inputfile, nb_div):
    line = None
    with open(inputfile, 'r') as fd:
        for line_num, line_content in enumerate(fd):
            if nb_div in line_content and line_num != j + 1:
                line = line_num
    return line


def FindtheCrispr(indexname, DR, count, seqbeg, seqend, crisprfile):
    spacers, goodcrispr, crisBeg, crisEnd, refFalsSpacers, simDRs = [], 0, 0, 0, [], None
    DRlength = len(DR)
    simDRs, refFalsSpacers, spacers = definespacers(DRlength, crisprfile, indexname)
    FalsSpacers = refFalsSpacers
    if len(FalsSpacers) > 0:
        FalsSpacers = [h + seqbeg for h in FalsSpacers]
    spacersH = trans_struct_hach(seqbeg, spacers)
    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}/log.txt', 'a', f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\t'
                                                 f'Check short CRISPRs arrays and check spacers alignment '
                                                 f'using Muscle...')
    o = len(spacers) - 1
    if minNbSpacers <= 0:
        minNbSpacers = 1
    if len(spacers) >= (minNbSpacers - 1):
        if len(spacers) == 0:
            goodcrispr = check_short_crispr(DR, "spacers" + indexname)
        else:
            goodcrispr = checkspacersAlignMuscle("spacers" + indexname)
        if goodcrispr:
            crisBeg = seqbeg + spacers[0].Pos1 - DRlength
            crisEnd = seqbeg + spacers[-1].Pos1 + spacers[-1].Length + DRlength - 1
    r = len(FalsSpacers) - 1
    return goodcrispr, simDRs, crisBeg, crisEnd, spacersH, len(spacers), FalsSpacers


def stdev(data):
    if len(data) == 1:
        return 0
    average = sum(data) / len(data)
    sqtotal = 0
    for i in data:
        sqtotal += (average - i) ** 2
    std = (sqtotal / (len(data) - 1)) ** 0.5
    return std


def create_recap(RefSeq, nbrcris, OneSpacerCris_nbr, ResultDir):
    directory = ResultDir
    Crispr_report = RefSeq + "_CRISPRs"

    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}/log.txt', 'a',
                   f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tCreate first CRISPR(s) file(s) ...........\n')

    if not os.path.isdir(directory):
        os.mkdir(directory)
    dir = directory + "/" + RefSeq
    if not os.path.isdir(dir):
        os.mkdir(dir)

    File_Content = "########################################\n"
    File_Content += "# Program: CRISPR Finder \n"
    File_Content += "# Author: Ibtissem GRISSA\n"
    File_Content += "# Rundate (GMT): %s\n" % date_time
    File_Content += "# Report file: %s\n" % Crispr_report
    File_Content += "########################################\n"
    File_Content += "#=======================================\n"
    File_Content += "# \n"
    File_Content += "# Sequence: %s \n" % RefSeq
    File_Content += "# Number of CRISPR-like structures : %s\n" % nbrcris
    File_Content += "# Number of questionable structures : %s\n" % OneSpacerCris_nbr
    File_Content += "# \n"
    File_Content += "#=======================================\n"
    os.chdir(dir)
    with open(Crispr_report, "w") as writer:
        writer.write(File_Content)


def trim(string):
    return string.strip()


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
               '-CASFinder': 'casfinder_opt', '-cf': 'casfinder_opt',
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
    bool_options = ['keep', 'logOption', 'html', 'kingdom', 'classifySmall', 'seqMinSize', 'mismOne', 'writeFullReport']

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
            if item == '-':
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
    if parametrs["repeatsQuery"] == '...':
        parametrs["repeatsQuery"] = ''
    num_options = ['SpSim', 'Sp2', 'Sp1', 'S2', 'S1', 'M2', 'M1', 'DRtrunMism', 'DRerrors', 'vicinity',
                   'flankingRegion', 'genCode', 'levelMin', 'percentageMismatchesHalfDR']
    for item in num_options:
        try:
            parametrs[item] = int(parametrs[item])
        except ValueError:
            try:
                parametrs[item] = float(parametrs[item])
            except ValueError:
                print(f'Error in {num_options}')

    # Проверка наличия входного файла
    if not os.path.isfile(parametrs['userfile']):
        sys.exit(f'Not found {parametrs["userfile"]}. '
                 f'Please check that the file exists or enter a correct file name.\n')

    # Коррекция DRs
    parametrs['DRtrunMism'] = 100 / float(parametrs['DRtrunMism'])
    parametrs['DRerrors'] = float(parametrs['DRerrors']) / 100

    # отметка времени при запуске процесса
    print(f'Launch Time: {datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}')

    # запуск работы с BIO | to_do
    if os.path.splitext(parametrs["userfile"])[-1] == '.fastq':
        SeqIO.convert(parametrs["userfile"], "fastq", f'{os.path.splitext(parametrs["userfile"])[-2]}.fasta', "fasta")
        parametrs["userfile"] = os.path.splitext(parametrs["userfile"])[-2] + '.fasta'
    SeqIO.read(parametrs["userfile"], "fasta")
    inputfileCount = 0

    # создание папки с итоговыми результатами
    outDir = f'{os.getcwd()}/{parametrs["outputDirName"]}'
    outDir_tsv = f'{outDir}/TSV'
    if not os.path.isdir(outDir):
        os.mkdir(outDir)
    if not os.path.isdir(outDir_tsv):
        os.mkdir(outDir_tsv)

    outDir_gff = f'{outDir}/GFF'
    if not os.path.isdir(outDir_gff):
        os.mkdir(outDir_gff)

    # для логов (ПОСМОТРЕТЬ LOGGING)
    if json.loads(parametrs['keep'].lower()):
        outDir_log = f'{outDir}/Temporary File'
    if not os.path.isdir(outDir_log):
        os.mkdir(outDir_log)
    if json.loads(parametrs['logOption'].lower()) and json.loads(parametrs['keep'].lower()):
        with open(f'{outDir_log}/log.txt', 'w') as logfile:
            lof_write = csv.writer(logfile, delimiter='\t')
            lof_write.writerow([parametrs["userfile"], f'SIZE: {os.path.getsize(parametrs["userfile"])} bytes'])
        logfile.close()

    # JSON file 'result'    373
    JSONRES = f'{outDir}/jsonResult.json'
    write_text(JSONRES, 'w', '{\n')
    write_text(JSONRES, 'a', f'Date:  {datetime.now().strftime("%d-%m-%Y | %H:%M:%S")} '
                             f'[dd-mm-yy | hh-mm-ss]\nVersion:  python\n')
    write_text(JSONRES, 'a', f'Command:  {config["Launch Function"]["Function"]}\n')
    write_text(JSONRES, 'a', 'Sequences:\n[{')

    if json.loads(parametrs['logOption'].lower()):
        write_text(f'{outDir_log}/log.txt', 'a',
                   f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tResults will be stored in {outDir}\n')

    if not json.loads(parametrs['quiet'].lower()):
        print(f'{datetime.now().strftime("%H:%M:%S")} ---> Results will be stored in {outDir}\n')

    allFoundCrisprs = 0
    allCrisprs = 0
    nbrAllCas = 0

    actualMetaOptionValue = parametrs['metagenome']

    # разобраться с этими принтами
    write_text(f'{outDir_tsv}/CRISPR-Cas_summary.tsv', 'w',
               'Sequence(s)\tCRISPR array(s)\tNb CRISPRs\tEvidence-levels\tCas cluster(s)\tNb Cas\t'
               'Cas Types/Subtypes\n')

    if json.loads(parametrs['clusteringThreshold'].lower()) and json.loads(parametrs['launchCasFinder'].lower()):
        write_text(f'{outDir_tsv}/CRISPR-Cas_clusters.tsv', 'w',
                   'Sequence\tCluster Id\tNb CRISPRs\tNb Cas\tStart\tEnd\tDescription\n')

    if json.loads(parametrs['classifySmall'].lower()):
        write_text(f'{outDir}/smallArraysReclassification.xlsx', 'w', '')
        write_text(f'{outDir}/smallArraysReclassification.xlsx', 'a',
                   'CRISPR_Id\tConsensus_Repeat\tFormer_Evidence-level\tNew_Evidence-level\n')
        # файл xls, не работает

    # подвязать pandas к созданию XLS
    # if json.loads(parametrs['classifySmall'].lower()):
    #    with open(f'{outDir}\\smallArraysReclassification.xls', 'wt') as smallArrays:
    #        smallArrays.write('text')
    #    smallArrays.close()

    seqIO = SeqIO.parse(parametrs['userfile'], 'fasta')
    # цикл 433-700 | работает с BIO
    for seq in seqIO:
        inputfileCount += 1
        inputfile = parametrs["userfile"]
        seqID1 = seq.id
        seqDesc = seq.description
        if seqDesc == '':
            seqDesc = 'Unknown'

        sequenceVersion = seqID1
        seqLength = len(seq)

        metagenome = actualMetaOptionValue

        if not metagenome and seqLength < 100000:
            metagenome = 1
        if seqLength >= int(parametrs['seqMinSize'].lower()):
            globalSeq = seq.seq
            globalAT = atpercent(globalSeq)

            # seqID1 = []
            # if "|" in seqID1:
            #     seqID1 = seqID1.split("|")
            #     seqID1 = seqID1.pop()
            #     seqID1 = []
            #     seqID1 = seqID1.split(".")
            #     seqID1 = seqID1[0]
            # else:
            #     if "." in seqID1:
            #         seqID1 = seqID1.split(".")
            #         seqID1 = seqID1[0]

            inputfile = f'{seqID1}.fna'
            Bio.SeqIO.write(seq, inputfile, 'fasta')

            nbrcris = 0
            OneSpacerCris_nbr = 0

            indexname = inputfile.split('.')
            print(f'Sequence number {inputfileCount}..\n')

            callmkvtree(inputfile, indexname)

            if not json.loads(parametrs['quiet'].lower()):
                print(f'( Input file: {inputfile}, Sequence ID: {seqID1}, Sequence name = {seqDesc} )\n')

            if os.path.isfile(parametrs['so']):
                if json.loads(parametrs['mismOne'].lower()):
                    if parametrs['repeatsQuery'] != '':   # os.path.isfile(parametrs['repeatsQuery'])
                        print('1')
                        vmatchoptions = f'-l {parametrs["M1"]} -s leftseq -evalue {1} -absolute -nodist' \
                                        f'-noevalue -noscore -noidentity -sort ida -q {parametrs["repeatsQuery"]}' \
                                        f'-best {1000000} -selfun {parametrs["so"]} {parametrs["M2"]}'
                    else:
                        print('2')
                        vmatchoptions = f'-l {parametrs["M1"]} {parametrs["S1"]} {parametrs["S2"]} -s leftseq -evalue' \
                                        f' {1} -absolute -nodist -noevalue -noscore -noidentity -sort ida -best' \
                                        f' {1000000} -selfun {parametrs["so"]} {parametrs["M2"]}'
                else:
                    print('3')
                    if parametrs['repeatsQuery'] != '':
                        vmatchoptions = f'-l {parametrs["M1"]} -e {parametrs["mismOne"]} -s leftseq -evalue {1}' \
                                        f' -absolute -nodist -noevalue -noscore -noidentity -sort ida -q' \
                                        f' {parametrs["repeatsQuery"]} -best {1000000} -selfun {parametrs["so"]}' \
                                        f' {parametrs["M2"]}'
                    else:
                        print('4')
                        vmatchoptions = f'-l {parametrs["M1"]} {parametrs["S1"]} {parametrs["S2"]} -e' \
                                        f' {parametrs["mismOne"]} -s leftseq -evalue {1} -absolute -nodist' \
                                        f' -noevalue -noscore -noidentity -sort ida -best {1000000} -selfun' \
                                        f' {parametrs["so"]}, {parametrs["M2"]}'
            else:
                print(f'The shared object file ({parametrs["so"]}) must be available in your current directory. '
                      f'Otherwise, you must use option -soFile (or -so)! \n\n')
                if json.loads(parametrs['logOption'].lower()):
                    write_text(f'{outDir_log}/log.txt', 'a',
                               f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tThe shared object file '
                               f'({parametrs["so"]}) must be available in your current directory. '
                               f'Otherwise, you must use option -soFile (or -so)! \n')
            vmatchoptions += f' {inputfile} > {outDir_tsv}/vmatch_result.txt'

            if json.loads(parametrs['logOption'].lower()):
                write_text(f'{outDir_log}/log.txt', 'a',
                           f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\tvmatch2 {vmatchoptions}\n')
            # =======================================================================================
            # =======================================================================================
            print('Контрольная точка открыта')
            # =======================================================================================
            makesystemcall(f'vmatch2 {vmatchoptions}')
            # =======================================================================================
            print('Контрольная точка закрыта')
            # =======================================================================================
            # =======================================================================================
            rep = trans_data(f'{outDir_tsv}/vmatch_result.txt')
            RefSeq = seqID1
            if len(rep) >= 0:
                nbrcris, OneSpacerCris_nbr = write_clusters(RefSeq, *rep)
                create_recap(RefSeq, nbrcris, OneSpacerCris_nbr, parametrs["outputDirName"])
                # outDir
            else:
                create_recap(RefSeq, 0, 0, parametrs["outputDirName"])

            actual_pathHome = getcwd()
            if json.loads(parametrs['logOption'].lower()):
                write_text(f'{outDir_log}/log.txt', 'a', f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\t '
                                                         f'Actual path in Home is: {actual_pathHome}\n')
            os.chdir("..")
            os.unlink("..\/crispr_result_*")
            os.unlink("..\/seq_v*")
            os.unlink(f"..\/{inputfile}.*")
            os.unlink("..\/spacersseq_v*")
            os.unlink("..\/fastaMuscle_spacersseq_v*")

            os.unlink('../stdout')
            os.unlink('../vmatch_result.txt')
            os.unlink('../alignDR_Spacer.needle')
            os.unlink('../DR')
            os.unlink('*_CRISPRs')
            os.chdir("..")
    #
    #         makesystemcall(f'rm -f {inputfile}.index')
    #
    #         gffFilename, *idDir = makeGff(ResultDir, inputfile, seqDesc, nbrcris, OneSpacerCris_nbr)
    #
    #         jsonFile = makeJson(gffFilename, ResultDir, RefSeq)
    #
    #         if not json.loads(parametrs['quiet'].lower()):
    #             print(f'Nb of CRISPRs in this sequence = {nbrcris}')
    #
    #         if json.loads(parametrs['logOption'].lower()):
    #             write_text(f'{outDir_log}/log.txt', 'a',
    #                        f'\n{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\t '
    #                        f'Nb of CRISPRs in this sequence {RefSeq} = {str(nbrcris)} \n')
    #
    #         analysis, foundCrisprs = crisprAnalysis(ResultDir, RefSeq, nbrcris, seqDesc, *idDir)
    #
    #         allFoundCrisprs += foundCrisprs
    #         allCrisprs += nbrcris
    #
    #         end_runCRISPR = datetime.now().strftime("%d-%m-%Y | %H:%M:%S")
    #
    #         casfile = ''
    #         jsonCAS = ''
    #         tabCRISPRCasClustersA = []
    #         nbrCas = 0
    #
    #         if rcfowce:
    #             if json.loads(parametrs['logOption'].lower()) and nbrcris > 0:
    #                 nbrCas, casfile, jsonCAS, *tabCRISPRCasClustersA = casFinder(ResultDir, inputfile, seqDesc, RefSeq,
    #                                                                              nbrcris, kingdom,
    #                                                                              parametrs['casfinder_opt'])
    #         else:
    #             if json.loads(parametrs['logOption'].lower()) and parametrs['rcfowce'] == 0:
    #                 nbrCas, casfile, jsonCAS, *tabCRISPRCasClustersA = casFinder(ResultDir, inputfile, seqDesc, RefSeq,
    #                                                                              nbrcris, kingdom,
    #                                                                              parametrs['casfinder_opt'])
    #
    #         if not nbrCas:
    #             nbrCas = 0
    #
    #         if not json.loads(parametrs['quiet'].lower()):
    #             print(f'Nb of Cas in this sequence = {nbrCas}')
    #
    #         if json.loads(parametrs['logOption'].lower()):
    #             write_text(f'{outDir_log}/log.txt', 'a',
    #                        f'{datetime.now().strftime("%d-%m-%Y | %H:%M:%S")}\t '
    #                        f'Nb of Cas in this sequence {RefSeq} = {str(nbrcris)} \n')
    #         nbrAllCas += nbrCas
    #
    #         # 625 - 631 | html
    #         # 631 - 702 | доделать
    #
    #         analyzedSequences = f'{outDir}/analyzedSequences'

    write_text(JSONRES, 'a', ']\n')
    write_text(JSONRES, 'a', '}\n')

    if json.loads(parametrs['launchCasFinder'].lower()) and json.loads(parametrs['writeFullReport'].lower()):
        fullReport(outDir)

    # if json.loads(parametrs['dirRepeat'].lower()):
    #     orientationCountFile = countOrientation(parametrs["outputDirName"], allCrisprs)
    #     if not json.loads(parametrs['quiet'].lower()):
    #         print("Orientations count file created: $orientationCountFile\n\n")
    # 740 - 773 | копия отчётов
    # 775 - 805 | результаты
    # 805 - 890 | конец программы


# подключение базовых настроек в INI файле
config = configparser.RawConfigParser()
config.optionxform = str
config.read("settings.ini")

parametrs = {}
[parametrs.update({item[0]: item[1]}) for item in config['Base Variable'].items()]

# запуск проверки наличия программ | to_do
print(f'Welcome to {config["System Variable"]["casfinder"]}.\n')
list_programm = ['vmatch2', 'mkvtree2', 'vsubseqselect2', 'fuzznuc', 'needle']
for item in list_programm:
    if isProgInstalled(item):
        print(f'_____{item} is found_____')
    else:
        sys.exit(f'Not found {item}')

active()  # основная ветка действий
