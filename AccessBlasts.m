function [ accession ] = AccessBlasts( NCBI,N )

RID = blastncbi(NCBI,'blastn');

disp('Waiting for reponse from NCBI blast')
pause(1);
disp('Wait for it...')
pause(1);
disp('Wait...')
pause(28); %if you do it too fast, the blast won't be up yet.

blast = getblast(RID);
cblast={blast.Hits.Name};
nblast=cblast(1:N)';

for ii=1:N;
sblast=strsplit(nblast{ii},'|');
access(ii)=sblast(4);
end
accession=access'
end

