import numpy as np
from Bio import SeqIO
import math
import pandas as pd
import csv
from keras.models import load_model

def countnum(seq,nuacid):
    return len([1 for i in range(len(seq)) if seq.startswith(nuacid,i)])

def construct_kmer():
	ntarr = ("A","C","G","T")

	kmerArray = []


	for n in range(4):
		kmerArray.append(ntarr[n])

	for n in range(4):
		str1 = ntarr[n]
		for m in range(4):
			str2 = str1 + ntarr[m]
			kmerArray.append(str2)
#############################################
	for n in range(4):
		str1 = ntarr[n]
		for m in range(4):
			str2 = str1 + ntarr[m]
			for x in range(4):
				str3 = str2 + ntarr[x]
				kmerArray.append(str3)
#############################################
#change this part for 3mer or 4mer
	for n in range(4):
		str1 = ntarr[n]
		for m in range(4):
			str2 = str1 + ntarr[m]
			for x in range(4):
				str3 = str2 + ntarr[x]
				for y in range(4):
					str4 = str3 + ntarr[y]
					kmerArray.append(str4)
############################################
	for n in range(4):
		str1 = ntarr[n]
		for m in range(4):
			str2 = str1 + ntarr[m]
			for x in range(4):
				str3 = str2 + ntarr[x]
				for y in range(4):
					str4 = str3 + ntarr[y]
					for z in range(4):
						str5 = str4 + ntarr[z]
						kmerArray.append(str5)
####################### 6-mer ##############
	for n in range(4):
		str1 = ntarr[n]
		for m in range(4):
			str2 = str1 + ntarr[m]
			for x in range(4):
				str3 = str2 + ntarr[x]
				for y in range(4):
					str4 = str3 + ntarr[y]
					for z in range(4):
						str5 = str4 + ntarr[z]
						for t in range(4):
							str6 = str5 + ntarr[t]
							kmerArray.append(str6)
####################### 7-mer ##############
	kmer7 = []
	for m in kmerArray[1364:5460]:
		for i in ntarr:
			st7 = m + i
			kmer7.append(st7)
	kmerArray = kmerArray+kmer7
    
	return kmerArray

def get_kmer(seq,kmerArray):
	rst = []
	total = 0.0
	for n in range(len(kmerArray)):
		item = countnum(seq,kmerArray[n])
		total = total + item
		rst.append(item) 
	for n in range(len(rst)):
		if total!=0:
			rst[n] = rst[n]/total

	return rst

def kmer_encode(seq,kmerArray):
    seq_data = []
    for i in range(len(seq)):
        seq_feature = get_kmer(seq[i],kmerArray)
        seq_data.append(seq_feature)
        if i==5000 or i==10000 or i==20000:
            print (i)
    seq_data = np.array(seq_data)
    return seq_data



# single nucleic ggap
def g_gap_single(seq,ggaparray,g):
    # seq length is fix =23
   
    rst = np.zeros((16))
    for i in range(len(seq)-1-g):
        str1 = seq[i]
        str2 = seq[i+1+g]
        idx = ggaparray.index(str1+str2)
        rst[idx] += 1
        
    for j in range(len(ggaparray)):
        rst[j] = rst[j]/(len(seq)-1-g) #l-1-g
        
    return rst

def ggap_encode(seq,ggaparray,g):
    result = []
    for x in seq:
        temp = g_gap_single(x,ggaparray,g)
        result.append(temp)
    result = np.array(result)
    return result

# binucleic ggap
# kmerarray[64:340]
def big_gap_single(seq,ggaparray,g):
    # seq length is fix =23
   
    rst = np.zeros((256))
    for i in range(len(seq)-1-g):
        str1 = seq[i]+seq[i+1]
        str2 = seq[i+g]+seq[i+1+g]
        idx = ggaparray.index(str1+str2)
        rst[idx] += 1
        
    for j in range(len(ggaparray)):
        rst[j] = rst[j]/(len(seq)-1-g) #l-1-g
        
    return rst

def biggap_encode(seq,ggaparray,g):
    result = []
    for x in seq:
        temp = big_gap_single(x,ggaparray,g)
        result.append(temp)
    result = np.array(result)
    return result

def orf_single(seq):
    startss = []
    stopss = []
    starts_ = []
    stop_s = []
    start = 0
    sslen = 0
    s_len1 = 0
    s_len2 = 0
    newseq = seq
    newseq_6 = []
    max_l = len(seq)
    l = len(seq)
    for i in range(len(seq)):
        if (seq[i:i+3]=="ATG"):    
            start = 1 # has start codon 
            for j in range(int((len(seq)-(i+3))/3)):
                if (seq[i+3+3*j:i+3+3*j+3]=="TAA") or (seq[i+3+3*j:i+3+3*j+3]=="TAG") or (seq[i+3+3*j:i+3+3*j+3]=="TGA"):
                    startss.append(i)
                    stopss.append(i+3+3*j+3)
                    break
            if len(startss)==0 :
                starts_.append(i)
                
    if start == 0:
        for k in range(len(seq)):
            if (seq[k:k+3]=="TAA") or (seq[k:k+3]=="TAG") or (seq[k:k+3]=="TGA"):
                stop_s.append(k+3)
        
    if len(startss)!=0:
        startss = np.array(startss)
        stopss = np.array(stopss)
        coding_len = stopss-startss
        max_len_position = np.argmax(coding_len)
        sslen = coding_len[max_len_position]
        newseq = seq[(startss[max_len_position]):(stopss[max_len_position])]
        max_l = sslen
        if (startss[max_len_position]-3)>=0 and (startss[max_len_position]+5)<l:
            newseq_6 = seq[(startss[max_len_position]-3):       (startss[max_len_position])]+seq[(startss[max_len_position]+3):(startss[max_len_position]+6)]                   
        
    elif len(starts_)!=0:
        starts_ = np.array(starts_)
        s_len1 = len(seq)-starts_[0]
        newseq = seq[(starts_[0]):len(seq)]
        max_l = s_len1
        if (starts_[0]-3)>=0 and (starts_[0]+5)<l:
            newseq_6 = seq[(starts_[0]-3):(starts_[0])]+seq[(starts_[0]+3):(starts_[0]+6)]           
              
    elif len(stop_s)!=0:
        stop_s = np.array(stop_s)
        s_len1 = stop_s[-1]
        newseq = seq[0:(stop_s[-1])]
        max_l = s_len1   
    
    orf_feature = (sslen/len(seq),s_len1/len(seq))

    return orf_feature,max_l,newseq,newseq_6

def orf_feature(seq):
    orf = []
    max_l = []
    newseq = []
    newseq_nu6 = []
    for i in range(len(seq)):
        orfsin,max_lsin,newseqsin,newseq_nu6sin = orf_single(seq[i])
        orf.append(orfsin)
        max_l.append(max_lsin)
        newseq.append(newseqsin)
        newseq_nu6.append(newseq_nu6sin)
    orf = np.array(orf)
    max_l = np.array(max_l)
    return orf,max_l,newseq,newseq_nu6

def fickett_single(seq):
    l = len(seq)
    a = []
    c = []
    g = []
    t = []
    for i in range(3):
        seq_3 = list(seq[i:l:3])
        a.append(seq_3.count('A'))
        c.append(seq_3.count('C'))
        g.append(seq_3.count('G'))
        t.append(seq_3.count('T'))
    a = np.array(a)    
    c = np.array(c)    
    g = np.array(g)    
    t = np.array(t)   
    apos = max(a)/(min(a)+1)
    cpos = max(c)/(min(c)+1)
    gpos = max(g)/(min(g)+1)
    tpos = max(t)/(min(t)+1)
    return (apos,cpos,gpos,tpos)        
def fickett_score(seq):
    rst = []
    for i in range(len(seq)):
        temp = fickett_single(seq[i])
        rst.append(temp)
    rst = np.array(rst)
    return rst

def hexamer_sin(seq,nc_m,c_m,kmerarray,x):
        
    if len(seq)>5:
        l = len(seq)-x+1
        log_r = np.zeros((l))
        for i in range(l):
            tempseq = seq[i:i+x]
            idx = kmerarray.index(tempseq)
            Fc = c_m[int(idx)]
            Fnc = nc_m[int(idx)]
            if Fc==0 and Fnc==0:
                log_r[i]=0
            elif Fc==0 and Fnc!=0:
                log_r[i]=-1
            elif Fnc==0 and Fc!=0:
                log_r[i]=1
            else:
                log_r[i] = math.log(Fc/Fnc)
        miu = sum(log_r)/l
    else:
        miu=0
    
    return miu
   
def hexamer_score(seq,nc_m,c_m,kmerarray,x):
    miu = np.zeros((len(seq)))
    for i in range(len(seq)):
        miu[i] = hexamer_sin(seq[i],nc_m,c_m,kmerarray,x)
        if i==5000 or i==10000 or i==20000:
            print (i)
    miu0 = np.expand_dims(miu, axis=1)
    return miu0

def nuc_count(seq):
    num=0
    for i in range(len(seq)):
        if len(seq[i])==6:
            num+=1
    nu5 = np.zeros((6,4))
    for i in range(6):
        temp = []
        for j in range(len(seq)):
            if len(seq[j])==6:
                temp.append(seq[j][i])
        nu5[i,:] = np.array([temp.count('A')/num,temp.count('C')/num,temp.count('G')/num,temp.count('T')/num])
    return nu5

def nucleic_biasin(seq,arrlnc,arrpc):
    if len(seq)==6:
        bia = np.zeros((6))
        for i in range(6):
            if seq[i]=='A':
                bia[i]=math.log(arrpc[i,0]/arrlnc[i,0])
            elif seq[i]=='C':     
                bia[i]=math.log(arrpc[i,1]/arrlnc[i,1])
            elif seq[i]=='G':      
                bia[i]=math.log(arrpc[i,2]/arrlnc[i,2])
            elif seq[i]=='T':        
                bia[i]=math.log(arrpc[i,3]/arrlnc[i,3])  
        bia0 = sum(bia)/6
    else:
        bia0=0
    return bia0
      
def nucleic_bia(seq,arrlnc,arrpc):
    bia=np.zeros((len(seq)))
    for i in range(len(seq)):
        bia[i] = nucleic_biasin(seq[i],arrlnc,arrpc)
    bia0 = np.expand_dims(bia, axis=1)
    return bia0

def encode(datapath,RNA_seq,m_type):
    if m_type=='human':
        max_fic=[11. , 10.6 ,12.5 ,13.9]
        maxl=107976
        lnc_arr=np.array([[0.28170942, 0.21526643, 0.26515386, 0.23787029],
                          [0.3228555 , 0.19672421, 0.27208512, 0.20833517],
                          [0.2809589 , 0.27989934, 0.28877312, 0.15036864],
                          [0.32016247, 0.20511236, 0.29945698, 0.1752682 ],
                          [0.28660986, 0.24943711, 0.22087325 ,0.24307978],
                          [0.24012185, 0.22952629, 0.28674231 ,0.24360955]])
        pc_arr= np.array([[0.43128673 ,0.12353085, 0.35409525, 0.09108717],
                          [0.29952253, 0.34491308, 0.21244491, 0.14311949],
                          [0.20277302, 0.42718536, 0.28121939, 0.08882223],
                          [0.23105411, 0.15548482, 0.46709721, 0.14636386],
                          [0.28281097, 0.37524486, 0.17911361, 0.16283056],
                          [0.17296156, 0.24978575, 0.35489104, 0.22236166]])

        nc_m = pd.read_csv(datapath+'humanlnc_6mermean.csv',header=None,delimiter = ',')
        c_m = pd.read_csv(datapath+'humanmrna_6mermean.csv',header=None,delimiter = ',')
        nc_m = np.array(nc_m)
        c_m = np.array(c_m)
    elif m_type=='vert':
        max_fic=[8.21978022 ,12.83333333 , 6.35 ,    7.85]
        maxl=100404
        lnc_arr=np.array([[0.27473965, 0.21007579, 0.25576577, 0.25941879],
                          [0.32642713, 0.20451448, 0.24595169, 0.2231067 ],
                          [0.26476201, 0.30238264, 0.28155499, 0.15130037],
                          [0.27664795, 0.24687858, 0.28215474, 0.19431874],
                          [0.29147811, 0.2575105 , 0.2074587 , 0.2435527 ],
                          [0.23924541, 0.21596423, 0.27626629, 0.26852407]])
        pc_arr= np.array([[0.51748992 ,0.08219269 ,0.33700324, 0.06331416],
                          [0.3221583 , 0.33144879, 0.18795874, 0.15843417],
                          [0.20921775, 0.42769292, 0.27653243, 0.0865569 ],
                          [0.22998082, 0.13800172, 0.48022879, 0.15178867],
                          [0.26942406, 0.41102956, 0.17020432, 0.14934206],
                          [0.15717781, 0.23999868, 0.36675924, 0.23606427]])
        nc_m = pd.read_csv(datapath+'vertlnc_6mermean.csv',header=None,delimiter = ',')
        c_m = pd.read_csv(datapath+'vertmrna_6mermean.csv',header=None,delimiter = ',')
        nc_m = np.array(nc_m)
        c_m = np.array(c_m)
    elif m_type=='insect':
        max_fic=[34.5 ,  19.82352941, 14.29032258,  9.81818182]
        maxl=48672
        lnc_arr=np.array([[0.30817248, 0.2036961 , 0.20320329, 0.28492813],
                          [0.32386037, 0.21084189, 0.24895277, 0.21634497],
                          [0.34373717, 0.21396304, 0.26743326, 0.17486653],
                          [0.31564682, 0.22776181, 0.24377823, 0.21281314],
                          [0.30784394, 0.22579055, 0.1949076 , 0.27145791],
                          [0.28082136, 0.21782341, 0.24320329, 0.25815195]])
        pc_arr= np.array([[0.61003182, 0.08203106, 0.23692585, 0.07101126],
                          [0.43908232, 0.27939011, 0.14102503, 0.14050254],
                          [0.39809053, 0.26461787, 0.228566  , 0.1087256 ],
                          [0.26499786, 0.17869187, 0.34398898, 0.21232128],
                          [0.29235738, 0.33358666, 0.17835938, 0.19569658],
                          [0.19080416, 0.2744977 , 0.29938726, 0.23531088]])
        nc_m = pd.read_csv(datapath+'insectlnc_6mermean.csv',header=None,delimiter = ',')
        c_m = pd.read_csv(datapath+'insectmrna_6mermean.csv',header=None,delimiter = ',')
        nc_m = np.array(nc_m)
        c_m = np.array(c_m)
    else:
        print ("type error")
    kmerArray = construct_kmer() 
    
    RNA_orf,RNA_maxlen0,RNA_seqnew,nuc6 = orf_feature(RNA_seq)

    RNA_hexamer = hexamer_score(RNA_seqnew,nc_m,c_m,kmerArray[1364:5460],6)

    RNA_fickett0 = fickett_score(RNA_seq)
    RNA_fickett = np.zeros((len(RNA_fickett0),4))

    RNA_fickett[:,0] = RNA_fickett0[:,0]/max_fic[0]
    RNA_fickett[:,1] = RNA_fickett0[:,1]/max_fic[1]
    RNA_fickett[:,2] = RNA_fickett0[:,2]/max_fic[2]
    RNA_fickett[:,3] = RNA_fickett0[:,3]/max_fic[3]
    
    RNA_maxlen1 = RNA_maxlen0/maxl
    RNA_maxlen = np.expand_dims(RNA_maxlen1, axis=1)

    RNA_biggap2 = biggap_encode(RNA_seq,kmerArray[84:340],2)
    RNA_biggap3 = biggap_encode(RNA_seq,kmerArray[84:340],3)

    RNA_kmer1 = kmer_encode(RNA_seq,kmerArray[0:4])
    RNA_kmer3 = kmer_encode(RNA_seqnew,kmerArray[20:84])

    nucbia = nucleic_bia(nuc6,lnc_arr,pc_arr)
    
    RNA_data = np.concatenate((RNA_hexamer,RNA_orf,RNA_maxlen,nucbia,RNA_fickett,RNA_kmer1,
                                  RNA_biggap2,RNA_biggap3,RNA_kmer3),axis=1) #    
    return RNA_data

def test_model(datapath,outpath,datafile,m_type,d_type):
    if m_type=='human':
        if d_type=='normal':
            model_name=['human1.h5','human2.h5','human3.h5']
            featureset = 'normal566.csv'
        elif d_type=='sorf':
            model_name=['hsorf1.h5','hsorf2.h5','hsorf3.h5']
            featureset = 'sorf450.csv'
        else:
            print ("type error")
    elif m_type=='vert':
        if d_type=='normal':
            model_name=['vert1.h5','vert2.h5','vert3.h5']
            featureset = 'normal566.csv'
        elif d_type=='sorf':
            model_name=['vsorf1.h5','vsorf2.h5','vsorf3.h5']
            featureset = 'sorf450.csv'
        else:
            print ("type error")
    elif m_type=='insect':
        if d_type=='normal':
            model_name=['insect1.h5','insect2.h5','insect3.h5']
            featureset = 'normal566.csv'
        elif d_type=='sorf':
            model_name=['isorf1.h5','isorf2.h5','isorf3.h5']
            featureset = 'sorf450.csv'
        else:
            print ("type error")
    else:
        print ("species error")
    records = list(SeqIO.parse(datapath+datafile, "fasta")) #  Human_orf_test_nc
    RNA_seq = []
    RNA_id=[]
    for i in range(len(records)):
        RNA_seq.append(records[i].seq)   
        RNA_id.append(records[i].id)
    maxposition=pd.read_csv(datapath+featureset,header=None,delimiter = ',') #testcpp_fea535 normal566 sorf450
    maxposition = np.array(maxposition)
    maxposition=[int(maxposition[i]) for i in range(len(maxposition))]
    RNA_data0 = encode(datapath,RNA_seq,m_type)
    RNA_data = RNA_data0[:,maxposition]
    testdata = np.expand_dims(RNA_data, axis=2) 
    
    bs=256
    probas = np.zeros((len(testdata),2))
    for filename in model_name:
        model = load_model(datapath+filename)   
        pre_rst = model.predict(testdata) # predict predict_classes
        probas = probas + pre_rst
        print ("finish model")
    l=3
    probas = probas/l
    # positive data right number
    pos_num = 0
    # negative data right number
    neg_num = 0
    cp=[]
    for i in range(len(probas)):
        if (probas[i,0]>0.5):
            cp.append('Noncoding')         
        else:
            cp.append('Coding')
    
    with open(outpath+'predict_results.csv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerows(zip(RNA_id,cp))
    
    return
        
