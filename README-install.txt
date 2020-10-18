FRESCO   FRES version 3.4

Source at https://github.com/I-Thompson/fresco 
Extract in new directory by command:
  git clone https://github.com/I-Thompson/fresco .

This directory contains four sub-directories: source, man, test and util.

The source/ directory contains : Fortran files *.f, makefile

The test/   directory contains : at least 6 test jobs xeta, lane20 & f19xfr,
                                             e80f49b, on2 & be11, etc
                                 their various outputs  SUN/*.out 
               (The input files were originally CRAY UNICOS jobs,
                hence the comments at the beginning.)

The man/   directory contains the instruction manual in latex:
                      fres32.tex, fres32.pdf:  latex source  & output
		 and fresco9.pdf: the Computer Physics Reports paper.

To compile FRESCO,
 
   1) Enter frxy/source, and then edit the makefile for your target machine,
	by set FC to your compiler, and FFLAGS to your compiler options
	
	The script 'mk' looks for which compilers are available,
	and selects the 'best' as $COMP by being later in the script.
	It then initiates compile in a subdirectory of $ARCH-$COMP,
	so multiple compilers and architecturs do not clash with each other.

   2) Edit aliases there,
      to set FRESCOLIB to point to directory for storing the binary

   3) Copy your aliases to ~/.fresco.aliases
      Execute .fresco.aliases e.g. in .cshrc  by including:	
        source ~/.fresco.aliases
      (This works in csh and tcsh, but not bash)

   4) Compile by:
        make    (compile in the source directory)
	 or
	mk      (compile in the source/$ARCH-$COMP directory)

   5) Install, copy 'fresco' and 'sfresco' to the FRESCOLIB by
        make install
	 or
	mk install

	Otherwise, copy them to some directory which is in your PATH list

   6) Clean up, with:	
        make clean
	 or
	mk clean
  
To run FRESCO,

   1) Enter test/ directory.

   2) The scripts include commands to construct temporary 'data' files. 
	These scripts are run by just saying  e.g.
       lane20.job

   3) To save the output in a file .e.g. `out', run the scripts by
       lane20.job > out &
         or simply
       run lane20.job 
         to use input file lane20.job and produce output file lane20.out.
 	('run' is a csh/tcsh alias in the 'aliases' file)

   4) If you have separate `data' or `in' files, the command is
       fresco < lane20.in > lane20.out

   5) To save any other output files from fresco, e.g. fort.16 for 
      cross sections, 
         touch lane20.xsecs
         ln -s lane20.xsecs fort.16
      before running fresco.
      The file fort.16 may have to be called for016.dat on some machines.

Please let me know if you have any questions or problems:
    
   I-Thompson@llnl.gov

Cheers, Ian Thompson
December, 2012
