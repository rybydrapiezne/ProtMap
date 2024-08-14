import os

import pygame
import sys
from pymol import cmd, finish_launching
from consurfDB import ConsurfDB
import requests
from bs4 import BeautifulSoup


def get_Class(listOfPDBID):
    dir = "classes"
    integere = 0
    for prot in listOfPDBID:
        try:
            response = requests.get(f'https://www.rcsb.org/annotations/{prot}')
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print(f'failed to load data for {prot}')
        fileName = os.path.join(dir, prot)
        f = open(fileName, "w")
        data = response.text
        temp = []

        alphas = data.find("All alpha")
        betas = data.find("All beta")
        alphaBeta = data.find("Alpha Beta")
        alphaBeta2 = data.find("Alpha and beta")
        if alphas != -1:  # TODO find a better solution for this!
            if alphas < betas != -1:
                if alphas < alphaBeta != -1:
                    if alphas < alphaBeta2 != -1:
                        temp.append("alpha")
                    else:
                        if alphaBeta2 != -1:
                            temp.append("alpha-beta")
                        else:
                            temp.append("alpha")
                else:
                    if alphaBeta != -1:
                        temp.append("alpha-beta")
                    else:
                        temp.append("alpha")
            else:
                if betas != -1:
                    if betas < alphaBeta != -1:
                        if betas < alphaBeta2 != -1:
                            temp.append("beta")
                        else:
                            if alphaBeta2 != -1:
                                temp.append("alpha-beta")
                            else:
                                temp.append("beta")
                    else:
                        if alphaBeta != -1:
                            temp.append("alpha-beta")
                        else:
                            temp.append("beta")
                else:
                    if alphas < alphaBeta != -1:
                        if alphas < alphaBeta2 != -1:
                            temp.append("alpha")
                        else:
                            if alphaBeta2 != -1:
                                temp.append("alpha-beta")
                            else:
                                temp.append("alpha")
                    else:
                        if alphaBeta != -1:
                            temp.append("alpha-beta")
                        else:
                            temp.append("alpha")
        elif betas != -1:
            if betas < alphaBeta != -1:
                if betas < alphaBeta2 != -1:
                    temp.append("beta")
                else:
                    if alphaBeta2 != -1:
                        temp.append("alpha-beta")
                    else:
                        temp.append("beta")
            else:
                if alphaBeta != -1:
                    temp.append("alpha-beta")
                else:
                    temp.append("beta")
        elif alphaBeta != -1:
            temp.append("alpha-beta")
        elif alphaBeta2 != -1:
            temp.append("alpha-beta")

        if "Membrane" in data:
            temp.append("membrane")

        if len(temp) != 0:
            f.write(temp[0])
        else:
            f.write("N/A")
            integere += 1
        f.close()
    print(integere)


def get_PDB_files(listOfPDBID):
    dir = "pdbs"

    for prot in listOfPDBID:
        try:
            response = requests.get(f'https://files.rcsb.org/view/{prot}.pdb')
            response.raise_for_status()
            fileName = os.path.join(dir, prot)
            f = open(fileName + ".pdb", "w")
            f.write(response.text)
            f.close()
        except requests.exceptions.RequestException as e:
            print(f'failed to load data for {prot}')


def get_PDB_list():
    temp = ""
    try:
        response = requests.get("https://pdb101.rcsb.org/motm/motm-image-download", verify=False)
        response.raise_for_status()
        data = response.text
        soup = BeautifulSoup(data, 'html.parser')
        # print(soup)
    except requests.exceptions.RequestException as e:
        print("failed to load pdb list: ", e)
        return None
    data = data.split("\n")
    tab = []
    temp = ""
    for line in data:
        temp = line.find(".tif")
        if temp != -1:
            tab.append(line[temp - 4:temp])
    tab = set(tab)
    last = []
    for i in tab:
        if i[0].isnumeric():
            if not i[1].isnumeric():
                last.append(i)

    # print(last)
    return last


def download_fastas(listOfPDBId):  # use only if no fasta files is present
    dir = "fastas"
    for prot in listOfPDBId:
        try:
            response = requests.get(f'https://www.rcsb.org/fasta/entry/{prot}')
            response.raise_for_status()
            fileName = os.path.join(dir, prot)
            f = open(fileName + ".fasta", "w")
            f.write(response.text)
            f.close()
        except requests.exceptions.RequestException as e:
            print(f'failed to load data for {prot}')


def get_from_consurf(proteinid):
    temp = ""
    try:
        response = requests.get(f'https://consurfdb.tau.ac.il/scripts/chain_selection.php?pdb_ID={proteinid}',
                                verify=False)

        response.raise_for_status()
        data = response.text
        soup = BeautifulSoup(data, 'html.parser')
        if soup.find('input', {'id': 'unique_chain'}) != None:
            unique_chain = soup.find('input', {'id': 'unique_chain'}).get('value')
            # pdb_chain = soup.find('input', {'id': 'PDB_chain'}).get('value')
            temp += unique_chain
        else:
            return None
        if temp == "":
            select_box = soup.find('select', {'name': 'PDB_chain_select'})
            option_values = [option['value'] for option in select_box.find_all('option') if option['value']]
            for i in option_values:
                if len(i) != 0:
                    return i[2:]
        return temp
    except requests.exceptions.RequestException as e:
        print("failed to load data: ", e)
        return None


def get_pdb(proteinid):
    temp = ""
    try:
        response = requests.get(f'https://consurfdb.tau.ac.il/DB/{proteinid}/{proteinid}_consurf_summary.txt',
                                verify=False)

        response.raise_for_status()
        data = response.text
        return data
    except requests.exceptions.RequestException as e:
        print("failed to load data: ", e)
        return None


def get_protein_sequence(pdb_file):
    #    finish_launching()
    #    cmd.load(pdb_file)
    sequence = {}


#    for chain in cmd.get_chains():
#        chain_seq = cmd.get_fastastr(f'chain {chain}')
#        sequence[chain] = chain_seq
#    return sequence
cp = []
# fileToTest = "7lq6.pdb"
listOfPDBId = get_PDB_list()
# for id in listOfPDBId:
#    cp.append(get_from_consurf(id))
# print(cp)
download_fastas(listOfPDBId)  # use only for gathering data
get_PDB_files(listOfPDBId)  # use only for gathering data
p = 0
get_Class(listOfPDBId)  # use only for gathering data

# for actId in cp:
#    if actId==None:
#        continue
#    pdb1 = get_pdb(actId)
#    if pdb1==None:
#        continue
#    #print(pdb1)
#    f=open(f"s/{actId}",'w')
#    f.write(pdb1)
directory = "s"
for file in os.listdir(directory):
    temporaryRes = []
    f = os.path.join(directory, file)
    with open(f, "r") as data:
        for line in data:
            line = line.split()
            temporaryRes.append(line)
        temporaryRes = temporaryRes[15:]
        temporaryRes = temporaryRes[0:len(temporaryRes) - 5]
    # print(temporaryRes)
# cmd.reinitialize()
# cmd.load("7lq6.pdb")

# Show the structure
# cmd.show("ribbon")
# cmd.bg_color("white")

# cmd.set('ray_opaque_background', 0)  # Make the background transparent
# cmd.zoom('all', -15)  # Zoom in by the specified factor

# Save an image
# cmd.png("output_image.png", width=300, height=300, dpi=600, ray=1)

# Quit PyMOL

# print(get_protein_sequence(fileToTest))
# fileName = "alpha.fas"
# currSequence = []
# integ = 0
# statsFile = "4D6TV_consurf_summary.txt"
# dataClean = []
# with open(statsFile) as data:
#    for line in data:
#        line = line.split()
#        if integ > 79:
#            dataClean.append(line)
#        integ += 1

# print(dataClean)
# with open(fileName) as sequence:
#    for line in sequence:
#        if line.startswith('>'):
#            continue  # Skip the header lines
#        else:
#            line = line.strip("\n")
#            currSequence.append(line)  # Add each character to the list

# pygame.init()
# inputRects = []
# inputLetters = []

# screen = pygame.display.set_mode([800, 600])
# clock = pygame.time.Clock()
# basefont = pygame.font.Font(None, 32)
# counter = 0
# currY = 100
# scrollx = 0

# Generate rects for each character
# for i in currSequence[:1]:
#    for char in i:
#        inputRects.append(pygame.Rect(5 + (30 * counter), currY, 25, 32))
#        inputLetters.append(char)
#        counter += 1

# integ = 0
# editableLetters = [0] * len(currSequence[0])
# for i in dataClean:
#    if float(i[3]) <= -0.9:
#        editableLetters[integ] = 1
#    integ += 1

# color_active = pygame.Color('lightskyblue3')
# color_passive = pygame.Color('chartreuse4')
# color_inactive = pygame.Color('red')
# active = False
# active_rect_index = -1

# while True:
#     for event in pygame.event.get():
#         if event.type == pygame.QUIT:
#             pygame.quit()
#             sys.exit()
#         if event.type == pygame.MOUSEWHEEL:
#             if event.y == 1:
#                 scrollx = min(scrollx + 30, 0)
#             elif event.y == -1:
#                 scrollx = max(scrollx - 30, -(len(currSequence[0]) * 30 - 800))  # Scroll left
#         if event.type == pygame.MOUSEBUTTONDOWN:
#             for i, rect in enumerate(inputRects):
#                 adjusted_rect = rect.move(scrollx, 0)
#                 if adjusted_rect.collidepoint(event.pos):
#                     if editableLetters[i] == 1:
#                         active = True
#                         active_rect_index = i
#                         break
#             else:
#                 active = False
#                 active_rect_index = -1
#         if event.type == pygame.KEYDOWN:
#             if event.key == pygame.K_LEFT:
#                 scrollx = min(scrollx + 30, 0)  # Scroll right
#             elif event.key == pygame.K_RIGHT:
#                 scrollx = max(scrollx - 30, -(len(currSequence[0]) * 30 - 800))  # Scroll left
#             elif active:
#                 if event.key == pygame.K_BACKSPACE:
#                     inputLetters[active_rect_index] = ''
#                 elif event.unicode.isprintable():
#                     inputLetters[active_rect_index] = event.unicode
#
#     screen.fill((255, 255, 255))
#
#     for i, input_rect in enumerate(inputRects):
#         adjusted_rect = input_rect.move(scrollx, 0)
#         if active and i == active_rect_index:
#             color = color_active
#         elif editableLetters[i] == 1:
#             color = color_passive
#         else:
#             color = color_inactive
#         pygame.draw.rect(screen, color, adjusted_rect, 2)
#
#         text_surface = basefont.render(inputLetters[i], True, (0, 0, 0))
#         screen.blit(text_surface, (adjusted_rect.x + 5, adjusted_rect.y + 5))
#
#     img = pygame.image.load("output_image.png").convert()
#     screen.blit(img, (100, 200))
#     pygame.display.flip()
#     clock.tick(60)
