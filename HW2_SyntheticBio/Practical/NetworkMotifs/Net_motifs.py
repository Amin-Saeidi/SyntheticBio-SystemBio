import random
import statistics
from scipy import stats

def data_preparation(path):
    data = []
    file = open(path, "r")
    for line in file:
        line = line.replace("\n","")
        data.append(line.split(" "))
    file.close()
    return data

def prepare_network(coliInterFullVec):
    network = dict()
    for element in coliInterFullVec:
        key = element[0]
        if key not in network:
            network[key] = []
        network[key].append([element[1],element[2]])
    return network

def generating_random_network(V, E):
    random_network = []
    for i in range(0,E):
        v_1 = random.randint(1,V)
        v_2 = random.randint(1,V)
        act_rep = random.randint(1,3)
        random_network.append([str(v_1),str(v_2),str(act_rep)])
    return random_network

def NOA_netwrok_generator(full_network):
    NOA_network = []
    for element in full_network:
        if element[0] != element[1]:
            NOA_network.append(element)
    return NOA_network
    
def Output_generator(V,E,number_of_networks,reality):
    output = dict()
    NARs = []
    PARs = []
    DPFLs = []
    DNFLs = []
    FAN_OUTs = []
    CASCADEs = []
    FFLs = []
    
    for i in range(0,number_of_networks):
    
        full_network = generating_random_network(V,E)
        network_NOA = NOA_netwrok_generator(full_network)
        full_network = prepare_network(full_network)
        network_NOA = prepare_network(network_NOA)
        NARs.append(len(finding_NAR(full_network)))
        PARs.append(len(finding_PAR(full_network)))
        DPFLs.append(len(finding_DPFL(network_NOA)))
        DNFLs.append(len(finding_DNFL(network_NOA)))
        FANOUTs = finding_FAN_outs(network_NOA)
        FAN_OUTs.append(len(FANOUTs))
        CASCADEs .append(len(finding_Cascades(network_NOA)))
        FFLs.append(len(finding_FFL(network_NOA,FANOUTs)))
    
    output['NARs'] = [statistics.mean(NARs),statistics.variance(NARs),stats.ttest_1samp(NARs, popmean=reality[0], alternative='two-sided').pvalue ]
    output['PARs'] = [statistics.mean(PARs),statistics.variance(PARs),stats.ttest_1samp(PARs, popmean=reality[1], alternative='two-sided').pvalue ]
    output['DPFLs'] = [statistics.mean(DPFLs),statistics.variance(DPFLs),stats.ttest_1samp(DPFLs, popmean=reality[2], alternative='two-sided').pvalue ]
    output['DNFLs'] = [statistics.mean(DNFLs),statistics.variance(DNFLs),stats.ttest_1samp(DNFLs, popmean=reality[3], alternative='two-sided').pvalue ]
    output['FAN_OUTs'] = [statistics.mean(FAN_OUTs),statistics.variance(FAN_OUTs),stats.ttest_1samp(FAN_OUTs, popmean=reality[4], alternative='two-sided').pvalue ]
    output['CASCADEs'] = [statistics.mean(CASCADEs),statistics.variance(CASCADEs),stats.ttest_1samp(CASCADEs, popmean=reality[5], alternative='two-sided').pvalue ]
    output['FFLs'] = [statistics.mean(FFLs),statistics.variance(FFLs),stats.ttest_1samp(FFLs, popmean=reality[6], alternative='two-sided').pvalue ]

    return (output)


# NODE HAS RELATION WITH ITSELF WITH THE NUMBER OF 2
def finding_NAR(network):
    NARs = []
    for key in network.keys():
        for TF in network[key]:
            if key == TF[0] and (TF[1] in ["2","3"]):
                NARs.append([key,TF[0],TF[1]]) 
    return NARs

# NODE HAS RELATION WITH ITSELF WITH THE NUMBER OF 1
def finding_PAR(network):
    PARs = []
    for key in network.keys():
        for TF in network[key]:
            if key == TF[0] and (TF[1] in ["1","3"]):
                PARs.append([key,TF[0],TF[1]]) 
    return PARs

def finding_DPFL(network):
    DPFL_pairs = []
    keys = network.keys()
    for key in keys:
        for TF in network[key]:
            key_1 = TF[0] 
            act_rep = TF[1]
            if act_rep in ["1","3"]:
                if key_1 in keys:
                    for value in network[key_1]:
                        if value[0] == key and (act_rep in ["1","3"]):
                            DPFL_pairs.append([key,key_1])   
    return DPFL_pairs

def finding_DNFL(network):
    DNFL_pairs = []
    keys = network.keys()
    for key in keys:
        for TF in network[key]:
            key_1 = TF[0] 
            act_rep = TF[1]
            if act_rep in ["2","3"]:
                if key_1 in keys:
                    for value in network[key_1]:
                        if value[0] == key and (act_rep in ["2","3"]):
                            DNFL_pairs.append([key,key_1])   
    return DNFL_pairs

def finding_FAN_outs(network):
    FAN_outs = dict()
    keys = network.keys()
    for key in keys:
        arr = []
        for TF in network[key]:
            act_rep = TF[1]
            if act_rep in ["1","3"]:
                arr.append(TF[0])
        if len(arr)>=2:        
            FAN_outs[key] = arr
    return FAN_outs

def finding_Cascades(network):
    Cascades_triples = []
    keys = network.keys()
    for key in keys:
        for TF in network[key]:
            key_1 = TF[0] 
            act_rep = TF[1]
            if act_rep in ["1","3"]:
                if key_1 in keys:
                    for value in network[key_1]:
                        if value[1]in ["1","3"]:
                            Cascades_triples.append([key,key_1,value[0]])   
    return Cascades_triples

def finding_FFL(network,Fanouts):
    FFL_triples = []
    keys = Fanouts.keys()
    for key in keys:
        for key_1 in Fanouts[key]:
            for key_2 in Fanouts[key]:
                if key_1 != key_2:
                    if key_1 in keys:
                        for value in network[key_1]:
                            if (value[1] in ["1","3"]) and value[0] == key_2:
                                FFL_triples.append([key,key_1,key_2])
    return FFL_triples

coliInterFullVec = data_preparation("coliInterFullVec.txt")
coliInterFullNames = data_preparation("coliInterFullNames.txt")
coliInterNoAutoRegVec = data_preparation("coliInterNoAutoRegVec.txt")

Number_of_Nodes = len(coliInterFullNames)
Number_of_Edges = len(coliInterFullVec)

network_A = prepare_network(coliInterFullVec)
network_NOA = prepare_network(coliInterNoAutoRegVec)

NARs = finding_NAR(network_A)
PARs = finding_PAR(network_A)
DPFLs = finding_DPFL(network_NOA)
DNFLs = finding_DNFL(network_NOA)
FAN_OUTs = finding_FAN_outs(network_NOA)
CASCADEs = finding_Cascades(network_NOA)
FFLs = finding_FFL(network_NOA,FAN_OUTs)

reality = [len(NARs),len(PARs),len(DPFLs),len(DNFLs),len(FAN_OUTs),len(CASCADEs),len(FFLs)]
sim = Output_generator(Number_of_Nodes,Number_of_Edges,1000,reality)

print(reality)
print(sim)