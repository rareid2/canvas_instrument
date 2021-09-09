Author : L. Madrias
Maintenace : J-Y Brochot
Date : 22/06/2006
Project : DEMETER
Version : 1.5


What is the function of the tool "rd_dmt_n1" ?
-------------------------------------------------------------------
This software reads n1 files (apid 1129 to 1144).
- It is possible to select a n1 file and the blocks to read,
- then these blocks are printed, packet per packet, in an ascii file or on the screen.


How to launch "rd_dmt_n1" ?
------------------------------------------
With the IDL virtual machine : 
	- copy "rd_dmt_n1.sav",
	- open the virtual machine
	- and select the file "rd_dmt_n1.sav"

With IDL :
	- copy "rd_dmt_n1.pro" and "rd_dmt_n1_sub.pro" in the same directory
	- open IDL
	- and enter "main" in the IDL command input line


How to use the tool "rd_dmt_n1" ?
--------------------------------------------------
1 - Select, clicking on the 1st button "Browse", the directory where the n1
file to read is located.
2 - Select the blocks to print.
3 - If you want the data to be printed in an ascii file, select "Yes", then 
choose, clicking on the 2nd button "Browse",  the directory where you want 
the output file to be saved. Click on "Compute" to begin the writting in the file.
The name of the ouptut file cannot be chosen by the user.
4 - If you want the data to be printed in the screen, click on "Display".
You can read packets by entering its number AND CLICKING ON "Display" or by 
directly using the arrows.
5 - Click on "Close" to close the window.


What are the outputs of the tool "rd_dmt_n1" ?
-------------------------------------------------------------------
Two outputs can be generated : 
1 - data are printed packet per packet on the screen,
2 - an ascii file can be written : 
its name is : "[name of the n1 file without the ".dat" extension]_lec.dat"


What are the versions of the tool "rd_dmt_n1" ?
-------------------------------------------------------------------
V1.0 : 04/03/2005
	- first version
V1.1 : 15/03/2005
	- correction of a bug (the ascii file now contains all the packets)
V1.2 : 04/04/2005
	- addition of the main procedure so rd_dmt_n1.sav can be rename
V1.3 : 25/05/2005
	- correction of a bug when reading ISL Survey file.
V1.4 : 08/06/2005
	- It's now possible to read data block4 for ISL!
V1.5 : 22/06/2006
        - Change of the labels of the fields 10,11 and 12 of ISL blocks 4.

Remark : 
-------------
For any information concerning n1 files, please read the document "Data Product 
Description", D. Lagoutte, J.Y. Brochot, M. Parrot, ref. DMT-SP-9-CM-6054_LPC)
