; Project : microsatellite Demeter
; Function : read Demeter N1 files
; Author : Line Madrias
; Source program : read_n1 from Jean-Yves Brochot
; Date : 25/05/2005
; Inputs :
; Outputs :
; Files needed : rd_dmt_n1_sub.pro
; Version : 1.3


; Modifications :
;   04/03/2005 : version 1.1 :
; 		bug corrected : the ascii file always begins by the 1st packet.
;   04/04/2005 : version 1.2 :
; 		addition of the main procedure so the saved file can be renamed
;   25/05/2005 : version 1.3 :
;                Bug corrected : when reading ISL Survey file, the
;                software hanged.
;   22/06/2006 : version 1.5 :
;                Change of the labels of the fields 10,11 and 12 of ISL
;                blocks 4.

; ##################################################################
; PRO rd_dmt_n1 : main procedure

; ##################################################################
; PRO rd_dmt_n1_event, event : management of the events of the main procedure




; Compilation of the other file(s)
@rd_dmt_n1_sub


; ################################################################################
; Events management
PRO rd_dmt_n1_event, event

	WIDGET_CONTROL, event.id, GET_UVALUE = uval
	CASE uval OF

		;-------------------------------------------------------------------------
		; Input file event
		'IBROWSE' : BEGIN
			WIDGET_CONTROL, event.top, GET_UVALUE = st, /NO_COPY
			CD, CURRENT = curdir
			IF (st.file_input_unit NE -1) THEN BEGIN
				IF ( (FSTAT(st.file_input_unit)).OPEN EQ 1 ) THEN BEGIN
					FREE_LUN, st.file_input_unit
					st.file_input = ''
					st.dir_input = ''
				ENDIF
			ENDIF

			; Filename by browse only
			filename = DIALOG_PICKFILE(DIALOG_PARENT = event.top, $
				/READ, $
				FILTER = 'DMT_N1_*', $
				TITLE = 'Select DMT_N1 file', $
				PATH = curdir, $
				GET_PATH = newdir)

			IF (filename NE '') THEN BEGIN
				PRINT, 'GOTEM'
				PRINT, st
				st.file_input = filename
				PRINT, st.file_input
				st.dir_input = newdir

				; Set the filename in the widget text
				WIDGET_CONTROL, st.id_text_fi_select, SET_VALUE = filename, SET_TEXT_SELECT = [STRLEN(filename)]
				WIDGET_CONTROL, st.id_base_param, /SENSITIVE
				WIDGET_CONTROL, st.id_base_fo, /SENSITIVE
				WIDGET_CONTROL, st.id_fo_choose, /SENSITIVE
				WIDGET_CONTROL, st.id_base_do, /SENSITIVE

				; Open the data at the first record
				GET_LUN, a
				st.file_input_unit = a
				OPENR, st.file_input_unit, st.file_input, /SWAP_IF_LITTLE_ENDIAN

				pos = STRPOS(st.file_input, "DMT_N1_", /REVERSE_SEARCH)
				st.apid = STRMID(st.file_input, pos+7, 4)

			ENDIF

		WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
		END


		;-------------------------------------------------------------------------
		; Packets to print
		'PARAM' : BEGIN
			WIDGET_CONTROL, event.top, GET_UVALUE = st, /NO_COPY

			; Gets the values modified by user
			WIDGET_CONTROL, st.id_param_choose, GET_VALUE = uparam
			st.val_param = LONG(uparam)
			WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
		END


		;-------------------------------------------------------------------------
		; Possibility to record
		'FO_YES_NO' : BEGIN
			WIDGET_CONTROL, event.top, GET_UVALUE = st, /NO_COPY
			WIDGET_CONTROL, st.id_fo_choose, GET_VALUE = v_yes_no
			st.val_yes_no = LONG(v_yes_no)
			IF (st.val_yes_no EQ 0L) THEN BEGIN
				WIDGET_CONTROL, st.id_button_fo_browse, /SENSITIVE
				WIDGET_CONTROL, st.id_text_fo_select, /SENSITIVE
			ENDIF ELSE BEGIN
				WIDGET_CONTROL, st.id_button_fo_browse, SENSITIVE = 0
				WIDGET_CONTROL, st.id_text_fo_select, SENSITIVE = 0
				WIDGET_CONTROL, st.id_button_fo_compute, SENSITIVE = 0
			ENDELSE

			WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
		END


		;-------------------------------------------------------------------------
		; Choice of the directory where the output file is recorded
		'OBROWSE' : BEGIN
			WIDGET_CONTROL, event.top, GET_UVALUE = st, /NO_COPY
			odir = DIALOG_PICKFILE(DIALOG_PARENT = event.top, $
				/DIRECTORY, $
				TITLE = 'Select directory', $
				PATH = st.dir_input)
			st.dir_output = odir

			; Set the directory in the widget text
			WIDGET_CONTROL, st.id_text_fo_select, SET_VALUE = odir, SET_TEXT_SELECT = [STRLEN(odir)]
			IF (odir NE '') THEN WIDGET_CONTROL, st.id_button_fo_compute, /SENSITIVE

	        WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
		END


		;-------------------------------------------------------------------------
		; Recording of the output file
		'OCOMPUTE' : BEGIN
			WIDGET_CONTROL, event.top, GET_UVALUE = st, /NO_COPY

			; Set the output filename depending on the input one
			CASE !VERSION.OS_FAMILY OF
				'Windows' : slash = '\'
				'unix' : slash = '/'
			ENDCASE
			n_pos = STRPOS(st.file_input, slash, /REVERSE_SEARCH)
			t = STRLEN(STRMID(st.file_input, n_pos+1))
			filen = st.dir_output + STRMID(st.file_input, n_pos+1, t-4) + '_lec' + '.dat'

			; Record the data in the output file
			GET_LUN, a
			st.file_output_unit = a
			OPENW, a, filen, /SWAP_IF_LITTLE_ENDIAN
			read_n1, st.apid, st.file_input_unit, st.file_output_unit, st.val_param
			PRINT, st.val_param
			FREE_LUN, st.file_output_unit

	        WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
		END


		;-------------------------------------------------------------------------
		; Choice of a number packet to read
		'PACKET_NB' :


		;-------------------------------------------------------------------------
		; Reading of the first packet
		'BEGINNING' : BEGIN
			WIDGET_CONTROL, event.top, GET_UVALUE = st, /NO_COPY

			st.val_packet = 0
			OPENW, temp_unit, 'temp_file.dat', /SWAP_IF_LITTLE_ENDIAN, /GET_LUN
			output_unit = temp_unit
			print_n1, st, output_unit
			FREE_LUN, output_unit
			a =WIDGET_BASE(/ROW)
			XDISPLAYFILE, 'temp_file.dat', WTEXT = w, GROUP = a
			WIDGET_CONTROL, w, GET_VALUE = str
			WIDGET_CONTROL, st.id_text_do_display, SET_VALUE = str
			WIDGET_CONTROL, a, /DESTROY
			WIDGET_CONTROL, st.id_text_do_select, SET_VALUE = STRING(st.val_packet)
	        FILE_DELETE, 'temp_file.dat'

	        WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
		END


		;-------------------------------------------------------------------------
		; Reading of the previous packet
		'BEFORE' : BEGIN
			WIDGET_CONTROL, event.top, GET_UVALUE = st, /NO_COPY

			st.val_packet = st.val_packet-1L > 0L
			OPENW, temp_unit, 'temp_file.dat', /SWAP_IF_LITTLE_ENDIAN, /GET_LUN
			output_unit = temp_unit
			print_n1, st, output_unit
			FREE_LUN, output_unit
			a =WIDGET_BASE(/ROW)
			XDISPLAYFILE, 'temp_file.dat', WTEXT = w, GROUP = a
			WIDGET_CONTROL, w, GET_VALUE = str
			WIDGET_CONTROL, st.id_text_do_display, SET_VALUE = str
			WIDGET_CONTROL, a, /DESTROY
			WIDGET_CONTROL, st.id_text_do_select, SET_VALUE = STRTRIM(STRING(st.val_packet), 2)
	        WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
	        FILE_DELETE, 'temp_file.dat'

	        WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
		END


		;-------------------------------------------------------------------------
		; Reading of the following packet
		'AFTER' : BEGIN
			WIDGET_CONTROL, event.top, GET_UVALUE = st, /NO_COPY

			st.val_packet += 1L

			OPENW, temp_unit, 'temp_file.dat', /SWAP_IF_LITTLE_ENDIAN, /GET_LUN
			output_unit = temp_unit
			print_n1, st, output_unit
			FREE_LUN, output_unit
			a =WIDGET_BASE(/ROW)
			XDISPLAYFILE, 'temp_file.dat', WTEXT = w, GROUP = a
			WIDGET_CONTROL, w, GET_VALUE = str
			WIDGET_CONTROL, st.id_text_do_display, SET_VALUE = str
			WIDGET_CONTROL, a, /DESTROY
			WIDGET_CONTROL, st.id_text_do_select, SET_VALUE = STRTRIM(STRING(st.val_packet), 2)
	        WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
	        FILE_DELETE, 'temp_file.dat'

	        WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
		END


		;-------------------------------------------------------------------------
		; Reading of the last packet
		'ENDING' : BEGIN
			WIDGET_CONTROL, event.top, GET_UVALUE = st, /NO_COPY

			st.val_packet = 100000L

			OPENW, temp_unit, 'temp_file.dat', /SWAP_IF_LITTLE_ENDIAN, /GET_LUN
			output_unit = temp_unit
			print_n1, st, output_unit
			FREE_LUN, output_unit
			a =WIDGET_BASE(/ROW)
			XDISPLAYFILE, 'temp_file.dat', WTEXT = w, GROUP = a
			WIDGET_CONTROL, w, GET_VALUE = str
			WIDGET_CONTROL, st.id_text_do_display, SET_VALUE = str
			WIDGET_CONTROL, a, /DESTROY
			WIDGET_CONTROL, st.id_text_do_select, SET_VALUE = STRTRIM(STRING(st.val_packet), 2)
	        WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
	        FILE_DELETE, 'temp_file.dat'

	        WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
		END


		;-------------------------------------------------------------------------
		; Reading of the chosen packet
		'DISPLAY' : BEGIN
			WIDGET_CONTROL, event.top, GET_UVALUE = st, /NO_COPY
			WIDGET_CONTROL, st.id_text_do_select, GET_VALUE = str_nb

			; Control the value
			nb = LONG(str_nb)
			nb = nb > 0L
			st.val_packet = nb

			OPENW, temp_unit, 'temp_file.dat', /SWAP_IF_LITTLE_ENDIAN, /GET_LUN
			output_unit = temp_unit
			print_n1, st, output_unit
			FREE_LUN, output_unit
			a =WIDGET_BASE(/ROW)
			XDISPLAYFILE, 'temp_file.dat', WTEXT = w, GROUP = a
			WIDGET_CONTROL, w, GET_VALUE = str
			WIDGET_CONTROL, st.id_text_do_display, SET_VALUE = str
			WIDGET_CONTROL, a, /DESTROY
			WIDGET_CONTROL, st.id_text_do_select, SET_VALUE = STRTRIM(STRING(st.val_packet), 2)
	        WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
	        FILE_DELETE, 'temp_file.dat'

	        WIDGET_CONTROL, event.top, SET_UVALUE = st, /NO_COPY
		END

		;-------------------------------------------------------------------------
		; Display of the data
		'TEXT_DISPLAY' :


		;-------------------------------------------------------------------------
		; Closing of all the widgets
		'CLOSE' : BEGIN
			WIDGET_CONTROL, event.top, GET_UVALUE = st, /NO_COPY

			IF (st.file_input_unit NE -1) THEN $
				IF ( (FSTAT(st.file_input_unit)).OPEN EQ 1 ) THEN FREE_LUN, st.file_input_unit
			IF (st.file_output_unit NE -1) THEN $
				IF ( (FSTAT(st.file_output_unit)).OPEN EQ 1 ) THEN FREE_LUN, st.file_output_unit
			CLOSE, /ALL
			WIDGET_CONTROL, st.id_base_main, /DESTROY
			WIDGET_CONTROL, DEFAULT_FONT = st.cur_font
		END

	ENDCASE


END
; End of the events management
; ################################################################################




; ################################################################################
; rd_dmt_n1 procedure
PRO rd_dmt_n1

	;-------------------------------------------------------------------------
	; Version
	version = '1.5'
	PRINT, 'RD_DMT_N1 software, version: ', version

	;-------------------------------------------------------------------------
	; Get current font and set new font (to be compatible on Windows and Unix)
	DEVICE, GET_CURRENT_FONT = cfont
	IF !D.NAME EQ 'WIN' THEN WIDGET_CONTROL, DEFAULT_FONT = 'helvetica*14'
;	IF !D.NAME EQ 'X' THEN WIDGET_CONTROL, DEFAULT_FONT = 'helvetica*11'
	IF !D.NAME EQ 'X' THEN WIDGET_CONTROL, DEFAULT_FONT = 'helvetica*10'


	;-------------------------------------------------------------------------
	; Main base
	Base_main = Widget_Base(GROUP_LEADER=wGroup, $
		UNAME = 'Base_main', $
		XOFFSET = 5, YOFFSET = 5, $
		SCR_XSIZE = 680, $
		SPACE = 5, $
		TITLE = 'READ DEMETER N1, V' + version, /COLUMN)


	;-------------------------------------------------------------------------
	; Input file selection
	Base_FI_Select = Widget_Base(Base_main, $
		UNAME = 'Base_FI_Select', $
		FRAME = 1, $
		XOFFSET = 7, YOFFSET = 14, $
		SCR_XSIZE = 640, SCR_YSIZE = 75, $
		SPACE = 3)

	Label_FI_Select = Widget_Label(Base_FI_Select, $
		UNAME = 'Label_FI_Select', $
		XOFFSET = 250, YOFFSET = 2 , $
		SCR_XSIZE = 158, $
		/ALIGN_CENTER, $
		VALUE='INPUT FILE SELECTION:')

	Button_FI_Browse = Widget_Button(Base_FI_Select, $
		UNAME = 'Button_FI_Browse', $
		XOFFSET = 10, YOFFSET = 28, $
		SCR_XSIZE = 56,  $
		/ALIGN_CENTER, $
		VALUE='Browse', $
		UVALUE='IBROWSE')

	Text_FI_Select = Widget_Text(Base_FI_Select, $
		UNAME = 'Text_FI_Select', $
		XOFFSET = 84, YOFFSET = 28, $
		SCR_XSIZE = 565, $
		VALUE = '    Select a DMT_<apid>_<startdate>_<enddate> file')


	;-------------------------------------------------------------------------
	; Block selection
	Base_Param = Widget_Base(Base_main, $
		UNAME = 'Base_Param_Select', $
		SENSITIVE = 0, $
		FRAME = 1, $
		XOFFSET = 7, YOFFSET = 99, $
		SCR_XSIZE = 640, SCR_YSIZE = 80, $
		SPACE = 3)

	params = [' Block 1', ' Block 2', ' Block 3', $
		'Beginning of Block 4', 'Data of Block 4']
	Param_Choose = CW_Bgroup(Base_Param, $
		params, $
		UVALUE = 'PARAM', $
		/ROW, $
		SPACE = 40, $
		XOFFSET = 35, YOFFSET = 2, $
		/NONEXCLUSIVE, $
		LABEL_TOP = 'BLOCKS TO PRINT:       ')


	;-------------------------------------------------------------------------
	; Output file selection
	Base_FO = Widget_Base(Base_main, $
		UNAME = 'Base_FO', $
		SENSITIVE = 0, $
		FRAME = 1, $
		XOFFSET = 7, YOFFSET = 99, $
		SCR_XSIZE = 640, SCR_YSIZE = 125, $
		SPACE = 3)

	FO_Choose = CW_Bgroup(Base_FO, $
		['Yes', 'No'] , $
		UVALUE = 'FO_YES_NO', $
		/ROW, $
		SPACE = 40, $
		XOFFSET = 150, YOFFSET = 2, $
		/EXCLUSIVE, $
		LABEL_LEFT = 'PRINT DATA IN AN ASCII FILE:                ')


	Button_FO_Browse = Widget_Button(Base_FO, $
		UNAME = 'Button_FO_Browse', $
		XOFFSET = 10, YOFFSET = 50, $
		SCR_XSIZE = 56,  $
		/ALIGN_CENTER, $
		VALUE='Browse', $
		UVALUE='OBROWSE', $
		SENSITIVE = 0)

	Text_FO_Select = Widget_Text(Base_FO, $
		UNAME = 'Text_FO_Select', $
		SENSITIVE = 0, $
		XOFFSET = 84, YOFFSET = 50, $
		SCR_XSIZE = 565, $
		VALUE = '    Select the directory for the output data file')

	Button_FO_Compute = Widget_Button(Base_FO, $
		UNAME = 'Button_FO_Compute', $
		XOFFSET = 300, YOFFSET = 90, $
		SCR_XSIZE = 56,  $
		SENSITIVE = 0, $
		/ALIGN_RIGHT, $
		VALUE='Compute', $
		UVALUE='OCOMPUTE')


	;-------------------------------------------------------------------------
	; Data per packet Display
	Base_DO = Widget_Base(Base_main, $
		UNAME = 'Base_DO', $
		FRAME = 1, $
		XOFFSET = 7, YOFFSET = 99, $
		SCR_XSIZE = 580, SCR_YSIZE = 580, $
		/SCROLL, $
		SPACE = 3, $
		SENSITIVE = 0)

	Label_DO = Widget_Label(Base_DO, $
		UNAME = 'Label_DO', $
		XOFFSET = 230, YOFFSET = 2 , $
		SCR_XSIZE = 200, $
		/ALIGN_CENTER, $
		VALUE='DATA DISPLAY PER PACKET:')

	Label_DO_Select = Widget_Label(Base_DO, $
		UNAME = 'Label_DO_Select', $
		XOFFSET = 15, YOFFSET = 30 , $
		SCR_XSIZE = 120, $
		VALUE='Packet number :      ')

	Text_DO_Select = Widget_Text(Base_DO, $
		UNAME = 'Text_DO_Select', $
		/EDITABLE, $
		XOFFSET = 150, YOFFSET = 28, $
		SCR_XSIZE = 50,  $
		/ALIGN_CENTER, $
		VALUE='0', $
		UVALUE='PACKET_NB')

	Button_Display = Widget_Button(Base_DO, $
		UNAME = 'Button_Display', $
		XOFFSET = 250, YOFFSET = 28, $
		SCR_XSIZE = 50,  $
		VALUE='Display', $
		UVALUE='DISPLAY')

	Button_DO_Beginning = Widget_Button(Base_DO, $
		UNAME = 'Button_DO_Beginning', $
		XOFFSET = 360, YOFFSET = 28, $
		SCR_XSIZE = 35,  $
		VALUE='l<<', $
		UVALUE='BEGINNING')

	Button_DO_Before = Widget_Button(Base_DO, $
		UNAME = 'Button_DO_Before', $
		XOFFSET = 420, YOFFSET = 28, $
		SCR_XSIZE = 35,  $
		VALUE='<', $
		UVALUE='BEFORE')

	Button_DO_After = Widget_Button(Base_DO, $
		UNAME = 'Button_DO_After', $
		XOFFSET = 480, YOFFSET = 28, $
		SCR_XSIZE = 35,  $
		VALUE='>', $
		UVALUE='AFTER')

	Button_DO_Ending = Widget_Button(Base_DO, $
		UNAME = 'Button_DO_Ending', $
		XOFFSET = 540, YOFFSET = 28, $
		SCR_XSIZE = 35,  $
		VALUE='>>l', $
		UVALUE='ENDING')

	Text_DO_Display = Widget_Text(Base_DO, $
		UNAME = 'Text_DO_Display', $
		XOFFSET = 7, YOFFSET = 60, $
		SCR_XSIZE = 620,  $
		SCR_YSIZE = 500,  $
		/SCROLL, $
		/ALIGN_CENTER, $
		VALUE='', $
		UVALUE='TEXT_DISPLAY')


	;-------------------------------------------------------------------------
	; Close the window and clean
	Button_Close = Widget_Button(Base_main, $
		VALUE = 'Close  ', $
		UVALUE = 'CLOSE', $
		/ALIGN_RIGHT)


	;-------------------------------------------------------------------------
	; Structure definition (all the parameters for management and computation)
	st = { $
		id_base_main : Base_main, $

		id_base_fi_select : Base_FI_Select, $
		id_button_fi_browse : Button_FI_browse, $
		id_text_fi_select : Text_FI_Select, $

		id_base_param : Base_param, $
		id_param_choose : Param_Choose, $

		id_base_fo : Base_FO, $
		id_fo_choose : FO_Choose, $
		id_button_fo_browse : Button_FO_Browse, $
		id_text_fo_select : Text_FO_Select, $
		id_button_fo_compute : Button_FO_Compute, $

		id_base_do : Base_DO, $
		id_text_do_select : Text_DO_Select, $
		id_button_do_beginning : Button_DO_Beginning, $
		id_button_do_before : Button_DO_Before, $
		id_button_do_ater : Button_DO_After, $
		id_button_do_ending : Button_DO_Ending, $
		id_text_do_display : Text_DO_Display, $

		file_input : '', $
		dir_input : '', $
		file_output : '', $
		dir_output : '', $

		file_input_unit : -1L, $
		file_output_unit : -1L, $

		val_yes_no : 0L, $
		val_param : LONARR(N_ELEMENTS(params)), $
		val_packet : 0L, $

		apid : '', $

		cur_font : cfont}

	WIDGET_CONTROL, Base_main, SET_UVALUE = st


	;-------------------------------------------------------------------------
	; Realization
	WIDGET_CONTROL, Base_main, /REALIZE
	XMANAGER, 'rd_dmt_n1', Base_main, /NO_BLOCK

END
; End of rd_dmt_n1 procedure
; ################################################################################



; ################################################################################
; Begin of the main procedure
PRO main
	rd_dmt_n1
END
; End of the main procedure
; ################################################################################
