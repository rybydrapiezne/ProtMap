from transformers import T5Tokenizer, T5EncoderModel
import torch
import re
import numpy
import matplotlib
import os
import h5py



device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

tokenizer = T5Tokenizer.from_pretrained('Rostlab/ProstT5', do_lower_case=False)

model = T5EncoderModel.from_pretrained("Rostlab/ProstT5").to(device)

model.full() if device == 'cpu' else model.half()

fileBetas = []
betaNames = []
fileAlphas = []
alphaNames = []

integ = 0
directory = 'fastas'
flag=False
for file in os.listdir(directory):
    temporaryRes = []
    fileBetas=[]
    f = os.path.join(directory, file)

    with open(f, "r") as data:
        for line in data:
            if integ == 0:
                line = line.strip("\n")
                integ = 1
                if 'DNA' in line:
                    flag=True
                    continue
                betaNames.append(line)

            elif integ == 1:
                if flag==False:
                    line = line.strip("\n")
                    integ = 0
                    if 'U' in line:
                        continue
                    fileBetas.append(line)
                    break
                else:
                    flag=False
                    integ=0

    # with open("beta.fas") as file:
    #     for line in file:
    #         if integ==0:
    #             line=line.strip("\n")
    #             betaNames.append(line)
    #             integ=1
    #         elif integ==1:
    #             line=line.strip("\n")
    #             fileBetas.append(line)
    #             integ=0

    # integ=0
    # with open("alpha.fas") as file:
    #     for line in file:
    #         if integ==0:
    #             line=line.strip("\n")
    #             alphaNames.append(line)
    #             integ=1
    #         elif integ==1:
    #             line=line.strip("\n")
    #             fileAlphas.append(line)
    #             integ=0
    sequence_examples = fileBetas

    sequence_examples = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequence_examples]
    print(sequence_examples)
    # The direction of the translation is indicated by two special tokens:
    # if you go from AAs to 3Di (or if you want to embed AAs), you need to prepend "<AA2fold>"
    # if you go from 3Di to AAs (or if you want to embed 3Di), you need to prepend "<fold2AA>"
    sequence_examples = ["<AA2fold>" + " " + sequence_examples[0] if sequence_examples[0].isupper() else "<fold2AA>" + " " + sequence_examples[0]
                         ]
    print(sequence_examples)
    ids = tokenizer.batch_encode_plus(sequence_examples,
                                      add_special_tokens=True,
                                      padding="longest",
                                      return_tensors='pt').to(device)
    # generate embeddings
    with torch.no_grad():
        embedding_rpr = model(
            ids.input_ids,
            attention_mask=ids.attention_mask
        )

    # extract residue embeddings for the first ([0,:]) sequence in the batch and remove padded & special tokens, incl. prefix ([0,1:8])
    emb_0 = embedding_rpr.last_hidden_state[0]  # shape (7 x 1024)

    emb_0 = emb_0.detach().cpu().numpy()

    h5py.File(file[0:len(file) - 6] + ".h5", emb_0)
    #print(emb_0)
