CC = "C:/ti/ccsv7/tools/compiler/ti-cgt-c6000_8.2.2/bin/cl6x.exe" -mv6600 -O1 -g -I"./codegen/lib/getHomographyMatrix" -I"C:/ti/ccsv7/tools/compiler/ti-cgt-c6000_8.2.2/include" -c --c99 -pc -DDSP
AR = "C:/ti/ccsv7/tools/compiler/ti-cgt-c6000_8.2.2/bin/ar6x.exe" -ru
SRCFILES = ./codegen/lib/getHomographyMatrix/*.c
OBJFILES = ./*.obj
libgetHomographyMatrix.lib:$(SRCFILES)
	$(CC) $(SRCFILES)
	$(AR) libgetHomographyMatrix.lib $(OBJFILES)


