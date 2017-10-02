function [ human, nonhuman ] = CloseMatch( NCBI )

N=50;
Accession=AccessBlasts(NCBI,N);

%for detecting the closest match to humans
hcounter=1;
h=0;
while h==0
    hcheck=getgenbank(Accession{hcounter});
    hspecie=hcheck.Source;
    h=strfind(hspecie, 'Homo sapiens');
    hempty=isempty(h);
        if hempty==1
            h=0;
        end
    hcounter=hcounter+1;
    if hcounter==N
        break
    end
end


%For finding the closest non human different specie. It can distinguish
%even between different strains 
%esta me va a dar la más cercana diferente, ósea, si el input es humano, el
%más cercano al humano. Si el input es ratón, el más cercano a ratón.
counter=1;
x=1;

while x==1
check1=getgenbank(Accession{counter});
counter=counter+1;
check2=getgenbank(Accession{counter});
specie1=check1.Source;
specie2=check2.Source;
 x=strfind(specie2,specie1);
    empty=isempty(x);
    if empty==1
        x=0;
    end
    
    x=strfind(specie2, 'Homo sapiens');
    empty=isempty(x);
    if empty==1
        x=0;
    end
    if x>1
       x=0;
       counter=counter+1;
    end
    if counter==N
        break
    end
end

nonhuman=[{check2.Accession}, {check2.Definition}]

if h>0
   human=[{hcheck.Accession}, {hcheck.Definition}]
end
if h==0
    disp('No human match found')
end


end

