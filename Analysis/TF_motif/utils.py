import random
def dna_1hot(seq, seq_len=None, n_uniform=False):
    """dna_1hot
    Args:
      seq:       nucleotide sequence.
      seq_len:   length to extend/trim sequences to.
      n_uniform: represent N's as 0.25, forcing float16,
                 rather than sampling.
    Returns:
      seq_code: length by nucleotides array representation.
    """
    if seq_len is None:
        seq_len = len(seq)
        seq_start = 0
    else:
        if seq_len <= len(seq):
            # trim the sequence
            seq_trim = (len(seq) - seq_len) // 2
            seq = seq[seq_trim : seq_trim + seq_len]
            seq_start = 0
        else:
            seq_start = (seq_len - len(seq)) // 2
    seq = seq.upper()

    # map nt's to a matrix len(seq)x4 of 0's and 1's.
    if n_uniform:
        seq_code = np.zeros((seq_len, 4), dtype="float16")
    else:
        seq_code = np.zeros((seq_len, 4), dtype="bool")

    for i in range(seq_len):
        if i >= seq_start and i - seq_start < len(seq):
            nt = seq[i - seq_start]
            if nt == "A":
                seq_code[i, 0] = 1
            elif nt == "C":
                seq_code[i, 1] = 1
            elif nt == "G":
                seq_code[i, 2] = 1
            elif nt == "T":
                seq_code[i, 3] = 1
            else:
                if n_uniform:
                    seq_code[i, :] = 0.25
                else:
                    ni = random.randint(0, 3)
                    seq_code[i, ni] = 1
    return seq_code

def pred_on_fasta(fa, model,batch_size):
    records = list(SeqIO.parse(fa, "fasta"))
    seqs = [str(i.seq) for i in records]
    seqs_1hot = np.array([dna_1hot(i) for i in seqs])
    X=seqs_1hot.astype(np.float32).swapaxes(-1,1)
    n_samples = X.shape[0]
    y_pred = []
    for i in tqdm(range(0, n_samples, batch_size)):#
        X_batch = X[i:i+batch_size, ...]
        X_batch = torch.from_numpy(X_batch).to(device)
        y_pred_batch = model.forward(X_batch).squeeze().cpu().data.numpy()
        y_pred_batch=y_pred_batch.astype(np.float32)
    
        ad1=anndata.AnnData(y_pred_batch)
        ad1.X=sparse.csr_matrix(ad1.X)
        y_pred.append(ad1)
    ad=anndata.concat(y_pred,join='outer')
    y_impute=ad.X.toarray()
    return y_impute