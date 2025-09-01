
   function Installmex_ssm(recompile)

   curdir = pwd;  
   fprintf(' current directory is:  %s\n',curdir);    
   if (nargin==0); recompile = 0; end
   computer_model = computer;
   matlabversion = sscanf(version,'%f');
   matlabversion = matlabversion(1);
   tmp = version('-release'); 
   matlabrelease = str2num(tmp(1:4));
%%    
   fsp = filesep;
   if strcmp(computer_model,'PCWIN')
      str0 = ['''',matlabroot,'\extern\lib\win32\lcc\'' '];
      if (exist(eval(str0),'dir')==7)
         str1 = ['''',matlabroot,'\extern\lib\win32\lcc\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win32\lcc\libmwblas.lib''  '];
      else
         str1 = ['''',matlabroot,'\extern\lib\win32\microsoft\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win32\microsoft\libmwblas.lib''  ']; 
      end
      libstr = [str1,str2];
   elseif strcmp(computer_model,'PCWIN64')
      str_lcc = ['''',matlabroot,'\extern\lib\win64\lcc\'' '];
      str_mingw64 = ['''',matlabroot,'\extern\lib\win64\mingw64\'' '];  
      str_microsoft = ['''',matlabroot,'\extern\lib\win64\microsoft\'' '];        
      if (exist(eval(str_lcc),'dir')==7)
         str1 = ['''',matlabroot,'\extern\lib\win64\lcc\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win64\lcc\libmwblas.lib''  '];
      elseif (exist(eval(str_mingw64),'dir')==7)
         str1 = ['''',matlabroot,'\extern\lib\win64\mingw64\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win64\mingw64\libmwblas.lib''  ']; 
      elseif (exist(eval(str_microsoft),'dir')==7)
         str1 = ['''',matlabroot,'\extern\lib\win64\microsoft\libmwlapack.lib''  '];
         str2 = ['''',matlabroot,'\extern\lib\win64\microsoft\libmwblas.lib''  '];
      else
         error(' no compiler found'); 
      end 
      libstr = [str1,str2];
   else
      libstr = '  -lmwlapack -lmwblas  '; 
   end
   if strfind(computer_model,'MAC')
      mexcmd = 'mex -largeArrayDims  -output ';
   else
      mexcmd = 'mex -O -largeArrayDims  -output ';
   end       
%%
%%
%%
   src = [curdir,fsp,'src',fsp,'mexfun']; 
   eval(['cd ','src',fsp,'mexfun']); 
   fprintf ('\n Now compiling the mexFunctions in:\n'); 
   fprintf (' %s\n',src);    

   fname{1} = 'mexbwsolve'; 
   fname{2} = 'mexCondat';
   fname{3} = 'mexDSYEVX';
   fname{4} = 'mexfwsolve';
   fname{5} = 'mexmat';
   fname{6} = 'mexMatvec';
   fname{7} = 'mexscale'; 
   fname{8} = 'mexsigma_update_Classic_Lasso_SSNAL';
   fname{9} = 'mexsmat';
   fname{10} = 'mexsvec';
   fname{11} = 'mextriang';
%%
   ext = mexext; 
   for k = 1:length(fname)
      existmex = exist([fname{k},'.',ext]); 
      if (existmex ~= 3) | (recompile)
         cmd([mexcmd,fname{k},'  ',fname{k},'.c  ',libstr]);  
      end
   end 
   cd .. 
   cd ..
%%***********************************************
   function cmd(s) 
   
   fprintf(' %s\n',s); 
   eval(s); 
%%***********************************************
