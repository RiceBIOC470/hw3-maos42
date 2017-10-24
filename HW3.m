% GB comments
1.	80 Need to write out the final alignment 
2a. 70 Be careful with the use of the function showalignment. Feeding your align into the function only gives a snippet of the entire possible coding sequence and therefore outputs a artificially high percent alignment. If you determine the length of ERK1.sequence, you will notice that its length is 1902 bps long and not 1506. 
2b. 70 Same issue as 2a. 
2c. 70 Same issue as 2a. 
3a 100 
3b. 100
3c. 100  	
Overall: 86


%HW3

%% Problem 1 - Smith-Waterman alignment
% Consider two sequences 'GTAATCC' and 'GTATCCG'

% Construct the scoring matrix for this with the parameters:
% match value = 2, mismatch value = -1, and gap penalty = -1. Use your
% solution to get the optimal alignment. If you prefer, it is acceptable to do this with
% pencil and paper, you can then take a snapshot of your solution and
% include it in your repository. 



%% Problem 2 - using the NCBI databases and sequence alignments

% Erk proteins are critical signal transducers of MAP kinase signaling.
% Accessions numbers for ERK1 (also called MAPK3) and ERK2 (also called MAPK1) human mRNA are NM_002746 and
% NM_002745, respectively. 

% Part 1. Perform an alignment of the coding DNA sequences of ERK1 and
% ERK2. What fraction of base pairs in ERK1 can align to ERK2? 

hERK1=getgenbank('NM_002746');
hERK2=getgenbank('NM_002745');

[score, align, start]=swalign(hERK1.Sequence,hERK2.Sequence,'Alphabet','nt','Showscore',true);

matches=count(align(2,:), '|');
total=length(align(2,:));
fraction=matches/total
showalignment(align);

%Miguel Angel: 1053/1506 = 69.99% of alignment

% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

phERK1=getgenpept(hERK1.CDS.protein_id);
phERK2=getgenpept(hERK2.CDS.protein_id);

[pscore, palign, pstart]=swalign(phERK1.Sequence,phERK2.Sequence);
pmatches=count(palign(2,:), '|');
ptotal=length(palign(2,:));
pfraction=pmatches/ptotal
showalignment(palign);


%Miguel Angel: 305/346 = 88.15% alignment

% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 

mERK1=getgenbank('X64605.1');
pmERK1=getgenpept(mERK1.CDS.protein_id);
mERK2=getgenbank('D10939.1');
pmERK2=getgenpept(mERK2.CDS.protein_id);

[hm1score, hm1align, hm1start]=swalign(hERK1.Sequence, mERK1.Sequence);
hm1matches=count(hm1align(2,:), '|');
hm1total=length(hm1align(2,:));
hm1fraction=hm1matches/hm1total
showalignment(hm1align);

[hm2score, hm2align, hm2start]=swalign(hERK2.Sequence, mERK2.Sequence);
hm2matches=count(hm2align(2,:), '|');
hm2total=length(hm2align(2,:));
hm2fraction=hm2matches/hm2total
showalignment(hm2align);

[phm1score, phm1align, phm1start]=swalign(phERK1.Sequence, pmERK1.Sequence);
phm1matches=count(phm1align(2,:), '|');
phm1total=length(phm1align(2,:));
phm1fraction=phm1matches/phm1total
showalignment(phm1align);

[phm2score, phm2align, phm2start]=swalign(phERK2.Sequence, pmERK2.Sequence);
showalignment(phm2align);
phm2matches=count(phm2align(2,:), '|');
phm2total=length(phm2align(2,:));
phm2fraction=phm2matches/phm2total
showalignment(phm2align);

%Miguel Angel:
%The ERK1 human vs mouse, the gene sequence is 89.74% similar.
%The ERK2 human vs mouse, the gene sequence is 75.22% simiar.

%The protein for ERK1 human vs mouse is 97.44% similar.
%The protein for ERK2 human vs mouse is 99.44% similar.

%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 

AccessBlasts('NM_002746',10);

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 

CloseMatch('NM_002746');

% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 

%Human gene: Homo sapiens leucine rich repeat kinase 2 (LRRK2), mRNA
%Human gene Accession: NM_198578.3
%Dengue gene: Dengue virus 4 isolate SG(EHI)D4/46870Y11 envelope protein gene, partial cds
%Dengue gene Accession: JN544417.1 

CloseMatch('NM_198578'); %for human gene
CloseMatch('JX024757'); %for foreign gene

%for the human gene:
%the closest non human: 'XM_003825726 Pan paniscus leucine-rich repeat kinase 2 (LRRK2)?' 
%the closest to the human genome: 'NM_198578 XM_058513 XM_930820 XM_930832 XM_936405 XM_942804' 
%
%the closest non human specie is the chimpanzee, which makes sense because
%is the closest specie to us humans, besides that gene is an enzyme
%involved in apoptosis in neuroblastoma cells for mainting a healthy neural population
%so that also accounts for the closeness to us. 
%As for the human blast, it is showing variances of the same gene which
%also makes sense, since I'm blasting a human gene against the human
%database the result should come quickly.

%for the foregin gene:
%'JX024757'    'Dengue virus 4 isolate EHI310A129CY10, complete genome.''
%No human match found.
%
%This results make sense as well. Since I am blasting a partial envelope
%protein gene of the 4th serotype of the dengue virus. The closest match
%will be the complete genome of the Dengue Virus serotype 4 in which the 
%partial envelope protein will be. And since it's a
%viral gene, it won't show up in the blast agains the human database.
%


