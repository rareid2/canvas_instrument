; i hate this language
; this works, need to add in a way to search thrugh the file treem and get some data
; then update python script to have the same file naming
; FINALLY do some statistics
; Compilation of the other file(s)
@rd_dmt_n1_sub
PRO rd_dmt_n1
	file_input_name = COMMAND_LINE_ARGS(COUNT=argc) 
	file_output_name = file_input_name + 'p' 

	PRINT, 'parsed', file_input_name

	; Open the data at the first record
	GET_LUN, a
	file_input_unit = a
	OPENR, file_input_unit, file_input_name, /SWAP_IF_LITTLE_ENDIAN

	GET_LUN, b
	file_output_unit = b
	OPENW, b, file_output_name, /SWAP_IF_LITTLE_ENDIAN
	read_n1, 1131, file_input_unit, file_output_unit, [1, 1, 1, 1, 1]
	
	FREE_LUN, file_output_unit
END