Papyrus README file
===================

last modified: Oct. 24, 2024
Papyrus version: 5.7
Total number of files: 14

===================

Papyrus is an aggregated dataset of small molecule bioactivities.

Changes since version 05.6
    - ChEMBL data was updated to ChEMBL version 34
    - data of the IUPHAR/BPS Guide to PHARMACOLOGY has been included
    - data of Pickett et al.'s publication on MMP-12 has been included (ACS Med Chem Lett. 2011 Jan 13; 2(1): 28–33. DOI: 10.1021/ml100191f)

The very high-quality Papyrus++ data is contained in the 
05.7++_combined_set_without_stereochemistry.tsv.xz file.

The other data is split into 3 tab-separated files:
    - 05.7_combined_set_without_stereochemistry.tsv.xz
        This file contains the highly standardized and normalized data.
        Although not containing stereochemistry, it provides the second
        highest quality of the data (after Papyrus++). This file can be readily used.
    - 05.7_combined_set_with_stereochemistry.tsv.xz
        This file contains standardized and normalized activity values
        but molecules are not standardized. It is of lower confidence than
        the data above although containing stereochemistry.
        By "molecules are not standardized" we mean that sources have
        applied different standardization processes and the reported
        stereochemistry may vary for the molecule the activities relate to.
        It is recommended to standardize the data before using it.
    - 05.7_combined_set_protein_targets.tsv.xz
        This file contains data about protein targets (e.g. sequence, organism,
        classification).

In addition to these, molecular structures are available in SDfile format:
    - 05.7_combined_2D_set_without_stereochemistry.sd.xz
        This file contains the unique molecules found in
        05.7_combined_set_without_stereochemistry.tsv.xz
    - 05.7_combined_3D_set_with_stereochemistry.sd.xz
        This file contains the unique molecules found in
        05.7_combined_set_with_stereochemistry.tsv.xz

Description of the data in Papyrus:
    The 05.7_combined_set_protein_targets.tsv.xz file contains the following columns:
        - target_id         A unique Papyrus protein identifier.
                            It results from the concatenation of
                            accessions and mutations
                            (e.g. P47747_WT or P10721_V559D_T670I).
        - TID               The original protein identifier.
                            It results from the concatenation of
                            the type of protein identifier and the
                            protein identifier used in the original dataset
                            (e.g. ChEMBL:CHEMBL2882 or HGNC:PKD1)
        - UniProtID         The UniProt identifier of the sequence.
        - Status            Status of the UniProt entry, either 'reviewed'
                            if entry comes from SwissProt or 'unreviewed'
                            if comes from TrEMBL.
        - Organism          Organism of the protein.
        - Classification    Protein classification as given by ChEMBL
                            (version 30). Levels are separated by "->".
                            Multiple classifications are separated by a
                            semicolon (";")
        - Length            Length of the protein sequence.
        - Sequence          Protein sequence including mutations.

    The 05.7_combined_set_without_stereochemistry.tsv.xz file contains:
        - Activity_ID           A unique Papyrus compound-protein identifier.
                                It results from the concatenation of the 
                                molecule InChI and protein unique identifier.
        - Quality               Quality assigned to the activity data
        - source                Sources of the activity values separated by
                                semicolons.
        - CID                   Compound identifiers (only original identifiers
                                for ExCAPE-DB and ChEMBL data).
        - SMILES                Standardized SMILES notation of the molecule.
        - connectivity          Connectivity of the molecule.
        - InChIKey              InChIKey of the molecule.
        - InChI                 InChI of the molecule.
        - InChI_AuxInfo         InChI_AuxInfo of the molecule.
        - target_id             Unique protein identifier, corresponds to the
                                "target_id" in 05.4_combined_set_protein_targets.tsv.xz.
        - TID                   The original protein identifier.
                                It results from the concatenation of
                                the type of protein identifier and the
                                protein identifier used in the original dataset
                                (e.g. ChEMBL:CHEMBL2882 or HGNC:PKD1)
        - accession             Protein UniProt Accession
        - Protein_Type          Either "WT" for wild-type sequences or a sequence
                                of the mutations underscore-separated.
        - AID                   Assay identifiers (only original identifiers
                                for ExCAPE-DB and ChEMBL data).
        - doc_id                First published or filed document characterizing the
                                compound-protein interaction. This does not always
                                correspond to the source of the reported activities.
        - Year                  Year the document referred to by "doc_id" was published
                                or filed.
        - all_docs_id           List of documents values in "pchembl_value" originate from.
        - all_years             List of years the documents in "all_docs_id" were filed or published.
        - type_IC50             Binary indicator of the types of the values in "pchembl_value"
                                semicolon-separated. If 0 then the activity is not an IC50,
                                if 1 then it is. If a unique value, then all activities have
                                the same type_IC50.
        - type_EC50             Binary indicator of the types of the values in "pchembl_value"
                                semicolon-separated. If 0 then the activity is not an EC50,
                                if 1 then it is. If a unique value, then all activities have
                                the same type_EC50.
        - type_KD               Binary indicator of the types of the values in "pchembl_value"
                                semicolon-separated. If 0 then the activity is not an KD,
                                if 1 then it is. If a unique value, then all activities have
                                the same type_KD.
        - type_Ki               Binary indicator of the types of the values in "pchembl_value"
                                semicolon-separated. If 0 then the activity is not an Ki,
                                if 1 then it is. If a unique value, then all activities have
                                the same type_Ki.
        - type_other            Binary indicator of the types of the values in "pchembl_value"
                                semicolon-separated. If 0 then the activity is IC50, EC50, Ki or KD,
                                if 1 then it is not. If a unique value, then all activities have
                                the same type_other.
        - Activity_class        Binary indication of activity if "pchembl_value" is empty.
        - relation              Comparison sign to indicate if "pchembl_value" is exact ("=")
                                or censored ("<", ">", "<=", ">=").
        - pchembl_value         Ensemble of log-transformed activity values available for then
                                compound-protein pair at the highest quality possible.
        - pchembl_value_Mean    Mean average of "pchembl_value" activities.
        - pchembl_value_StdDev  Standard deviation of "pchembl_value" activities.
        - pchembl_value_SEM     Standard error of the mean of "pchembl_value" activities.
        - pchembl_value_N       Number of "pchembl_value" activities for that compound-protein pair.
        - pchembl_value_Median  Median of "pchembl_value" activities.
        - pchembl_value_MAD     Mean absolute error of "pchembl_value" activities.

    The 05.7_combined_set_with_stereochemistry.tsv.xz file contains the same columns as the
    05.7_combined_set_without_stereochemistry.tsv.xz file with few exceptions:
        - connectivity is not present
        - SMILES is not standardized but contains stereochemistry.
    
    The 05.7_combined_2D_set_without_stereochemistry.sd.xz file contains
    connection tables of molecules with 2D coordinates and few properties:
        - SMILES        Standardized SMILES notation of the molecule in the
                        related tab-separated file.
        - connectivity  Connectivity of the molecule that should be used to
                        join this file with its related tab-separated file.
        - InChIKey      InChIKey of the molecule.
        - NEMAKey       NEMAKey  of the molecule (internally used to identify
                        molecules with same connectivity but differing tautomeric
                        state).

    The 05.7_combined_3D_set_with_stereochemistry.sd.xz file contains
    connection tables of molecules with 3D coordinates and the same properties
    as the 05.7_combined_2D_set_without_stereochemistry.sd.xz file with a few exceptions.
        - connectivity  is no present. InChIKey should be used in order to join 
                        this file to 05.4_combined_set_with_stereochemistry.tsv.xz.
        - Embedding     The embedding engine that was used to generate 3D coordinates.

Description of the molecular descriptors available in Papyrus:

    The 05.7_combined_2D_moldescs_CDDDs.tsv.xz file contains Continuous Data-Driven 
    Descriptors (CDDD; Winter et al. (2019) Chem. Sci. 10(6), 1692-1701, DOI: 10.1039/c8sc04175j)
    calculated from the 2D molecular structures in 05.7_combined_2D_set_without_stereochemistry.sd.xz.

    The 05.7_combined_2D_moldescs_ECFP6.tsv.xz file contains the RDKit extended connectivity
    Morgan fingerprints with radius 3 and 2048 bits calculated from the 2D molecular structures
    in 05.7_combined_2D_set_without_stereochemistry.sd.xz.

    The 05.7_combined_2D_moldescs_mold2.tsv.xz file contains the Mold2 descriptors
    (Hong et al. (2008). J. Chem. Inf. Model. 48(7), 1337-1344, DOI: 10.1021/ci800038f)
    calculated from the 2D molecular structures in 05.7_combined_2D_set_without_stereochemistry.sd.xz.

    The 05.7_combined_2D_moldescs_mordred2D.tsv.xz file contains Mordred descriptors
    (Moriwaki et al. (2018) J. Cheminf., DOI: 10.1186/s13321-018-0258-y) calculated from
    the 2D molecular structures in 05.7_combined_2D_set_without_stereochemistry.sd.xz.

    The 05.7_combined_3D_moldescs_E3FP.tsv.xz file contains extended three-dimensional fingerprint
    (E3FP; Axen et al. (2017) J. Med. Chem. 60(17), 7393-7409, DOI: 10.1021/acs.jmedchem.7b00696)
    calculated from 3D molecular structures in 05.7_combined_3D_set_with_stereochemistry.xd.xz.

    The 05.7_combined_3D_moldescs_mordred3D.tsv.xz file contains Mordred descriptors calculated
    from 3D molecular structures in 05.7_combined_3D_set_with_stereochemistry.xd.xz.

Description of the protein descriptors available in Papyrus:

    The 05.7_combined_protdescs_ProDEC.tsv.xz file contains a collection of protein 
    descriptors like various z-scales, VHSE and BLOSUM. Transformed using domain averages 
    across 50 domains and auto cross-covariances and physicochemical distance both with a 
    lag of 20 amino acids calculated from the protein sequences in 
    05.7_combined_set_protein_targets.tsv.xz.

    The 05.7_combined_protdescs_unirep.tsv.xz file contains UniRep protein descriptors
    (Alley et al. (2019) Nat. Methods 16(12), 1315-1322, DOI: 10.1038/s41592-019-0598-1)
    calculated from the protein sequences in 05.7_combined_set_protein_targets.tsv.xz.

Description of the metadata for Papyrus:

    The 05.7_additional_files.zip contains 3 files:
        - LICENSE.txt containing the license
        - data_size.json containing the number of entries for each of the datasets and descriptor sets 
          in Papyrus
        - data_types.json containing the data types of all columns occurring in the previously described sets
