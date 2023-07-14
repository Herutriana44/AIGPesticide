import random
import pickle
from rdkit import Chem
import pandas as pd
from rdkit.Chem import Descriptors
import random
from rdkit import Chem
from rdkit.Chem import rdDistGeom
import pandas as pd
import datetime

# Load the trained model from pickle file
with open('AIGPesticideDiscriminator.pkl', 'rb') as file:
    nb_classifier = pickle.load(file)

# Load the vectorizer from pickle file
with open('vectorizer.pkl', 'rb') as file:
    vectorizer = pickle.load(file)

class GAN:
    def augment_data(smiles, n_augmentations):
        """
        Menghasilkan n_augmentations SMILES tambahan dari suatu SMILES yang diberikan.
        """
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return None
        augmentations = []
        for i in range(n_augmentations):
            # Clone the original molecule
            new_m = Chem.Mol(m)
            
            # Generate a new conformation for the molecule
            rdDistGeom.EmbedMolecule(new_m, maxAttempts=10)
            
            # Generate a new SMILES from the modified molecule
            new_smiles = Chem.MolToSmiles(new_m)
            
            # Append the new SMILES to the list of augmentations
            augmentations.append(new_smiles)
            
        return augmentations

    def LD50(smiles):
        mol = Chem.MolFromSmiles(smiles)
        logP = Descriptors.MolLogP(mol)
        LD50 = (1 - logP) * 1000
        return LD50

    data = pd.read_csv('ActiveCompoundPersticide.csv')

    def generate_random_smiles(n):
        data = pd.read_csv('ActiveCompoundPersticide.csv')
        smiles_random = data['Substance'][random.randint(0, len(data)-1)]

        # lakukan augmentasi data
        smiles_augmented = GAN.augment_data(smiles_random, n)
        return smiles_augmented
        

    def save_smiles(smiles):
        with open('generated_smiles.txt', 'a') as file:
            file.write(smiles+'\n')

    def fix_smiles(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol)

    def main(num_molecules=10):
        result = pd.DataFrame(columns=['Substance', 'LD50','Time'])
        res = []
        for i in range(num_molecules):
            while True:
                smiles = GAN.generate_random_smiles(1)
                smiles = smiles[0]
                # predict using the trained model
                pred = nb_classifier.predict(vectorizer.transform([smiles]))
                if pred == 1:
                    print(smiles)
                    cek = Chem.MolFromSmiles(smiles)
                    if cek is None:
                        print(f"SMILES {smiles} tidak valid")
                        smiles = GAN.fix_smiles(smiles)
                        print(f"SMILES {smiles} sudah diperbaiki")
                    waktu = datetime.datetime.now()
                    result = result.append({'Substance': smiles, 'LD50': GAN.LD50(smiles), 'Time': waktu}, ignore_index=True)
                    res.append(smiles)
                    print(GAN.LD50(smiles))
                    break
        data = pd.read_csv('static/data/CompoundPersticide.csv')
        Res = pd.concat([data, result])
        # ubah kolom Time menjadi tipe datetime
        Res['Time'] = pd.to_datetime(Res['Time'])
        # urutan dari kolom Time berdasarkan waktu terbaru
        Res = Res.sort_values(by=['Time'], ascending=False)
        Res.to_csv('static/data/CompoundPersticide.csv', index=False)
        return result
