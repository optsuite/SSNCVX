function Installmex_ssm 

   computer_model = computer;
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);
%%
   if strcmp(computer_model,'PCWIN')
      str1 = ['''',matlabroot,'\extern\lib\win32\mingw64\libmwlapack.lib''  ']; 
      str2 = ['''',matlabroot,'\extern\lib\win32\mingw64\libmwblas.lib''  '];
      libstr = [str1,str2];     
   elseif strcmp(computer_model,'PCWIN64')
      str1 = ['''',matlabroot,'\extern\lib\win64\mingw64\libmwlapack.lib''  ']; 
      str2 = ['''',matlabroot,'\extern\lib\win64\mingw64\libmwblas.lib''  '];
      libstr = [str1,str2];  
   elseif strcmp(computer_model,'MACI64') || strcmp(computer_model,'MACA64')
      setenv('DYLD_LIBRARY_PATH', '/opt/homebrew/opt/lapack/lib')
      libstr = '-lmwlapack -lmwblas  '; 
   else
      libstr = '-lmwlapack -lmwblas  '; 
   end
   mexcmd = 'mex -O  -largeArrayDims  -output ';    

   fsp = filesep;

   curdir = pwd;  
   fprintf(' current directory is:  %s\n',curdir); 
%%
%% generate mex files in mexfun
%%
   clear fname


   fprintf ('\n Now compiling the mexFunctions in:\n'); 
   %%
   fname{1} = 'mexbwsolve';
   fname{2} = 'mexfwsolve';
   fname{3} = 'mexmat';
   fname{4} = 'mexsmat'; 
   fname{5} = 'mexsvec'; 
   fname{6} = 'mextriang'; 
   fname{7} = 'mexScalarxMatrix'; 
 
   for k = 1:length(fname)
       cmd([mexcmd,fname{k},'  ',fname{k},'.c  ',libstr]);  
   end      
   
   cmd(['mex -O -lmwlapack -lmwblas mexDSYEVX.c']);
   if strcmp(computer_model,'PCWIN')
      cmd([mexcmd,'mexeig','  ','mexeigwin','.c  ',libstr]);
   else
      cmd([mexcmd,'mexeig','  ','mexeig','.c  ',libstr]);
   end
   fprintf ('\n Compilation of mexFunctions completed.\n'); 
   cd .. 
   cd ..
   cd example
%%***********************************************
   function cmd(s) 
   
   fprintf(' %s\n',s); 
   eval(s); 
%%***********************************************