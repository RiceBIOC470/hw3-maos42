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
showalignment(align);

%Miguel Angel: 1053/1506 = 70% of alignment

% Part2. Perform an alignment of the aminoacid sequences of ERK1 and ERK2.
% What fraction of amino acids align?

phERK1=getgenpept(hERK1.CDS.protein_id);
phERK2=getgenpept(hERK2.CDS.protein_id);

[pscore, palign, pstart]=swalign(phERK1.Sequence,phERK2.Sequence);
showalignment(palign);

%Miguel Angel: 305/346 = 88% alignment

% Part 3.  Use the NCBI tools to get mRNA sequences for the mouse genes ERK1 and
% ERK2 and align both the coding DNA sequences and protein sequences to the
% human versions. How similar are they? 

mERK1=getgenbank('X64605.1');
pmERK1=getgenpept(mERK1.CDS.protein_id);
mERK2=getgenbank('D10939.1');
pmERK2=getgenpept(mERK2.CDS.protein_id);

[hm1score, hm1align, hm1start]=swalign(hERK1.Sequence, mERK1.Sequence);
showalignment(hm1align);
[hm2score, hm2align, hm2start]=swalign(hERK2.Sequence, mERK2.Sequence);
showalignment(hm2align);


[phm1score, phm1align, phm1start]=swalign(phERK1.Sequence, pmERK1.Sequence);
showalignment(phm1align);
[phm2score, phm2align, phm2start]=swalign(phERK2.Sequence, pmERK2.Sequence);
showalignment(phm2align);

%Miguel Angel:
%The ERK1 human vs mouse, the gene sequence is 90% similar.
%The ERK2 human vs mouse, the gene sequence is 75% simiar.

%The protein for ERK1 human vs mouse is 97% similar.
%The protein for ERK2 human vs mouse is 99% similar.

%% Problem 3: using blast tools programatically

% Part 1. Write a function that takes an NCBI accession number and a number N as input and
% returns a cell array of the accession numbers for the top N blast hits. 

%Miguel �ngel: I put the example from the previous problem, I also put a
%pause operation within the function, so it will take about 30 seconds to
%get the list of blast hits, because if you do not pause it, there is a
%chance that the blast isn't ready in the NCBI database and will throw
%back an error. 

AccessBlasts('NM_002746',5);

% Part 2. Write a function that takes an accession number as input, calls your function 
% from part 1, and returns two outputs - the closest scoring match in human DNA/RNA and the
% closest non-human match. Hint: see the "Source" or "SourceOrganism" field in the data
% returned by getgenbank. Make sure your function does something sensible
% if nothing human is found. 




% Part 3. Choose at least one gene from the human genome and one gene from
% another organism and run your code from part 2 on the relevant accession
% numbers. Comment on the results. 


