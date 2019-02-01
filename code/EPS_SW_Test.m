clear
clc

%miu=1
source_miu1_orig = [2	3	5	5	6	6	6	8	7	9	8	9	13	16	11	16	16	18	18	3	4	6	7	10	16	18	22	38	24	24	28	31	31	29	30	30	31	31	31	38	36	36	38	35	38	40	41	42	43	38	44	44	45	46	48	50	48	52	52	54	52	53	38	15	49	15	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	3	3	4	4	22	22	24	24	8	8	33	33	50	50	7	7	25	25	52	52	46	46	47	47	54	54	55	55	41	41	44	44	45	45	56	56	12	12	13	13	53	53	35	35	36	36	51	51	5	5	21	21	23	23	1	1	2	2	17	17	20	20	10	10	11	11	27	27	29	29	30	30	32	32	6	6	9	9	26	26	28	28	31	31	34	34	37	37	38	38	48	48	49	49	39	39	40	40	42	42	43	43	16	16	18	18	19	19	3	10];
target_miu1_orig = [1	2	3	4	5	10	12	7	10	8	10	11	12	13	16	18	19	17	21	24	22	31	28	36	48	52	23	22	23	25	25	26	27	28	29	36	34	35	37	31	32	33	35	36	37	38	38	38	38	44	45	46	46	47	47	47	51	50	53	52	55	56	14	14	15	44	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	115	126	115	126	115	126	115	126	119	131	119	131	119	131	118	127	118	129	118	129	140	141	140	141	140	141	140	141	138	139	138	139	138	139	138	139	121	122	121	122	121	122	133	134	133	134	133	134	116	125	116	125	116	125	113	114	113	114	113	114	113	114	120	129	120	129	120	129	129	130	129	130	129	130	117	128	117	128	117	128	128	132	128	132	128	132	135	142	135	142	135	142	135	142	136	137	136	137	136	137	136	137	123	124	123	124	123	124	143	144];
filename = 'Matrix_Gen_New.xlsx';
data = xlsread(filename,1);
[rows, columns] = size(data);
n1 = 56;
n2 = 88;

G_random = graph(source_miu1_orig, target_miu1_orig);
M = numedges(G_random); %Calculation of number of edges in the graph
L_random = log(n1+n2)/log(2*M/(n1+n2)); %Calculation of random network Characteristic Path Length
Gamma_random = 2*M/((n1+n2)^2); %Calculation of random network Characteristic Path Length

source_miu1 = source_miu1_orig;
target_miu1 = target_miu1_orig;

for col = 1 : 1 : 1000
    x = data(1:232, col);
    for i = 1:232
        if x(i,1) < 1 && x(i,1) > 0
            x(i,1) = 1;
        elseif x(i,1) < -1
            s_i = source_miu1(1,i);
            t_i = target_miu1(1,i);
            source_miu1(1,i) = t_i;
            target_miu1(1,i) = s_i;
            x(i,1) = -x(i,1);
        elseif x(i,1) > -1 && x(i,1) < 0
            s_j = source_miu1(1,i);
            t_j = target_miu1(1,i);
            source_miu1(1,i) = t_j;
            target_miu1(1,i) = s_j;
            x(i,1) = 1;
        end
    end
end

G_miu1 = graph(source_miu1, target_miu1, 1./x);
source_miu1 = source_miu1_orig;
target_miu1 = target_miu1_orig;

nn = numnodes(G_miu1);
[s,t] = findedge(G_miu1);
A = sparse(s,t,G_miu1.Edges.Weight,nn,nn);
A = (A+A') - eye(size(A,1)).*diag(A);
ccfs(:,col) = clustering_coefficients(A);

ccfs_miu1(:,col) = mean(ccfs(:,col)); %Calculation of Clustering Coefficient



d_miu1 = distances(G_miu1);
[L_miu1,efficiency,ecc,radius,diameter] = charpath(d_miu1);%Calculation of Characteristic Path Length

% FIRST SPLIT

adjA_drop = A;

for numnode = 1:nn
    adjA_drop(numnode,:) = 0;
    adjA_drop(:,numnode) = 0;
    % Currently find the following line caused a problem!
    %ccfs_drop = clustering_coefficients(adjA_drop);%Calculation of Clustering Coefficient after node removal
    %ccfs_drop_mean(1,col) = mean(ccfs_drop);
    G_miu1_drop = graph(adjA_drop);
    d_miu1_drop = distances(G_miu1_drop);
    %[L_miu1_drop,efficiency,ecc,radius,diameter] = charpath(d_miu1_drop); %Calculation of Characteristic Path Length after node removal
    %Delta_L(numnode,col) = (L_miu1-L_miu1_drop)/L_miu1; %Calculation of Characteristic Path Length drop
    %SWNS_drop_miu1(numnode,col) = (ccfs_drop_mean(1,col)/Gamma_random)/(L_miu1_drop/L_random); %Calculation of Small World-ness after node removal
end
%SWNS_miu1(:,col) = (ccfs_mesh(1,col)/Gamma_random)/(L_miu1/L_random); %Calculation of Small World-ness
%SWNS_drop_tot_miu1(:,col) = SWNS_drop_miu1; %Calculation of total Small World-ness after node removal
