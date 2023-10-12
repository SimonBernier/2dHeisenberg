clearvars
Ly = 4; Lx = [25 35 45]; v = [1.5:0.1:2 3:12];
h = 6; tau = [0.5 1 2 4];
maxDim = [512 800]; truncerr = 1E-8;

A=length(Lx); B = length(v); C=length(tau); D=length(maxDim);

for a=0:A-1
    for b=0:B-1
        for c=0:C-1
            for d=0:D-1
                
                runNumber = b*C*D + c*D + d;

                fileID = fopen(sprintf('input_2dHeis_mf_Ly_%d_Lx_%d_runN_%d',Ly,Lx(a+1),runNumber),'w');
                fprintf(fileID, 'input\n\t{\n');
                fprintf(fileID,sprintf('\tLy = %d\n',Ly));
                fprintf(fileID,sprintf('\tLx = %d\n',Lx(a+1) ));
                fprintf(fileID,sprintf('\tv = %0.2f\n',v(b+1)));
                fprintf(fileID,sprintf('\th = %0.2f\n',h));
                fprintf(fileID,sprintf('\ttau = %0.2f\n',tau(c+1)));
                fprintf(fileID,sprintf('\tmaxDim = %d\n',maxDim(d+1)));
                fprintf(fileID,sprintf('\ttruncE = %0.0e\n',truncerr));
                fprintf(fileID,'\tGSETDVP = 1\n');
                fprintf(fileID, '\t}\n');
                fclose(fileID);
            
            end
        end
    end
end

%%
Ly = 4; Lx = [25 35 45]; h = 6; tau = [2 4 6 round((Lx-1)/pi)];
maxDim = [512 800]; truncerr = 1E-8;

A=length(Lx); B = length(tau); C=length(maxDim);

for a=0:A-1
    for b=0:B-1
        for c=0:C-1
                
            runNumber = b*C + c;

            fileID = fopen(sprintf('input_2dHeis_uni_Ly_%d_Lx_%d_runN_%d',Ly,Lx(a+1),runNumber),'w');
            fprintf(fileID, 'input\n\t{\n');
            fprintf(fileID,sprintf('\tLy = %d\n',Ly));
            fprintf(fileID,sprintf('\tLx = %d\n',Lx(a+1) ));
            fprintf(fileID,sprintf('\th = %0.2f\n',h));
            fprintf(fileID,sprintf('\ttau = %0.2f\n',tau(b+1)));
            fprintf(fileID,sprintf('\tmaxDim = %d\n',maxDim(c+1)));
            fprintf(fileID,sprintf('\ttruncE = %0.0e\n',truncerr));
            fprintf(fileID,'\tGSETDVP = 1\n');
            fprintf(fileID, '\t}\n');
            fclose(fileID);
            
        end
    end
end

