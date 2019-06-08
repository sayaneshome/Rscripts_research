## Fetch PDB files and split to chain A only PDB files
ids <- c("3prg", "2vv3", "2q61","2q59","2q6s","2q6r","2q5s","2pob","2q5p","2p4y","2i4z","2i4p","2i4j","2g0h","2hfp","2g0g","2fvj","1zeo","1zgy","1wm0","1k74","1rdt","1fm9","1fm6")

files <- get.pdb(ids, split = TRUE, path = tempdir())
