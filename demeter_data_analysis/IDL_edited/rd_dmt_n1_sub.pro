; Project : microsatellite Demeter
; Function : read Demeter N1 files
; Author : Jean-Yves Brochot
; Date : 07/02/2005
; Main File : rd_dmt_n1.pro
; Version : 1.1


; ##################################################################
; FUNCTION TimeToS, t : convert time to seconds
; FUNCTION SToTime, tps : convert seconds to time

; ##################################################################
; PRO read_block1, unit, t0, output_unit, dh
; PRO read_block2, unit, output_unit, dh
; PRO read_block3, unit, output_unit, dh
; PRO read_block4_ubf_wave, unit, output_unit, dh, d
; PRO read_block4_ebf_wave, file_input_unit, output_unit, dh,d
; PRO read_block4_tbf_wave, unit, output_unit, dh, d
; PRO read_block4_spect, unit, output_unit, dh, d
; PRO read_block4_hf_wave, unit, output_unit, dh, d
; PRO read_block4_iap, unit, output_unit, dh, d
; PRO read_block4_isl, unit, te, output_unit, dh, d
; PRO read_block4_idp_burst, unit, output_unit, dh, d
; PRO read_block4_idp_survey, unit, TT, CC1, CC2, CC3, output_unit, dh, d
; PRO read_block4_rnf, unit, output_unit, dh, d
; 	all these procedures are called by read_n1 or print_n1

; ##################################################################
; PRO read_n1, apid, unit, output_unit, params : main procedure,
; 	called by rd_dmt_n1 (in rd_dmt_n1.pro)

; ##################################################################
; PRO print_n1, st, output_unit : main procedure, called by rd_dmt_n1




;------------------------------------------------------------------------------
; Function : TimeToS
; Description : Convert a time [yyyy, mm, dd, hh, mn, ss, mls] or an
;               array of time into number of seconds elapsed since a
;               fixed date.
; INPUTS   : t : an INTARR(7) or INTARR(7, N)
; OUTPUTS  : DOUBLE or DBLARR(N)
; KEYWORDS : none
; COMMONS  : none
;------------------------------------------------------------------------------

FUNCTION TimeToS, t

	tps = JULDAY(t(1, *), t(2,*), t(0,*)) - 2440588L ; Nb jours depuis 1/1/1970
	tps = tps * 86400.D + $
		t(3,*) * 3600.D + $
		t(4, *) * 60.D + $
		t(5, *) $
		+ t(6,*) / 1000.

	tps = REFORM(tps, N_ELEMENTS(TPS), /OVERWRITE)
	IF (N_ELEMENTS(tps) EQ 1) THEN tps = tps[0]
	RETURN, tps

END
; End of TimeToS
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of SToTime
FUNCTION SToTime, tps

	jj = LONG(tps / 86400.)
	CALDAT, jj + 2440588L, mo, da, year ; 2440588 = Nb jours depuis 1/1/1970
	reste = LONG((tps - jj * 86400.D) * 1000. + 0.5)
	hh = reste / 3600000UL
	reste = reste mod 3600000L
	mn = reste / 60000L
	reste = reste mod 60000L
	ss = reste / 1000L
	mls = reste mod 1000L
	t = FIX([[year], [mo], [da], [hh], [mn], [ss], [mls]])

	IF (SIZE(tps, /N_DIMENSIONS) GT 0) THEN $
		t = REFORM(t, [size(jj, /DIMENSIONS), 7])

	t = TRANSPOSE(t)

	RETURN, t
END
; End of SToTime
;------------------------------------------------------------------------------



;------------------------------------------------------------------------------
; Begin of read_block1
PRO read_block1, unit, t0, output_unit, dh

	IF (dh EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block1"
		PRINTF, output_unit, "======"
		P = 0B
		o1 = 0UL
		READU, unit, o1
		P = BYTE(ISHFT(o1, -24))
		PRINTF, output_unit, P, format='("P = 0x", Z2.2)'
		nbd = o1 AND 'FFFFFF'X
		PRINTF, output_unit, "nbd", nbd
		mls = 0UL
		READU, unit, mls
		PRINTF, output_unit, "mls = ", mls
		time = INTARR(7)
		READU, unit, time
		PRINTF, output_unit, "time", time
		t0 = TimeToS(time)
		orb = 0U
		READU, unit, orb
		PRINTF, output_unit, "orb = ", orb
		sorb = 0U
		READU, unit, sorb
		PRINTF, output_unit, "sorb = ", sorb
		sta = BYTARR(8)
		READU, unit, sta
		sta = STRING(sta)
		PRINTF, output_unit, "station = ", sta
		softv = 0UB
		softsv = 0UB
		READU, unit, softv
		READU, unit, softsv
		PRINTF, output_unit, "Soft version = ", STRING(FORMAT='(I1.1, ".", I1.1)', $
			softv, softsv)
		calv = 0UB
		calsv = 0UB
		READU, unit, calv
		READU, unit, calsv
		PRINTF, output_unit, "Cal version = ", STRING(FORMAT='(I1.1, ".", I1.1)', $
			calv, calsv)
	ENDIF ELSE BEGIN
		trash = BYTARR(38)
		READU, unit, trash
	ENDELSE
END
; End of SToTime
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block2
PRO read_block2, unit, output_unit, dh

	IF (dh EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block2"
		PRINTF, output_unit, "======"
		; field 1
		geocentric_lat = 0.
		READU, unit, geocentric_lat
		PRINTF, output_unit, "Geocentric_lat = ", geocentric_lat
		; field 2
		geocentric_long = 0.
		READU, unit, geocentric_long
		PRINTF, output_unit, "Geocentric_long = ", geocentric_long
		; field 3
		altitude = 0.
		READU, unit, altitude
		PRINTF, output_unit, "Altitude = ", altitude
		; field 4
		local_time = 0.
		READU, unit, local_time
		PRINTF, output_unit, "Local_Time = ", local_time
		;field 5
		geomagnetic_lat = 0.
		READU, unit, geomagnetic_lat
		PRINTF, output_unit, "Geomagnetic_lat = ", geomagnetic_lat
		;field 6
		geomagnetic_long = 0.
		READU, unit, geomagnetic_long
		PRINTF, output_unit, "Geomagnetic_long = ", geomagnetic_long
		;field 7
		Mlt = 0.
		READU, unit, Mlt
		PRINTF, output_unit, "Mlt = ", Mlt
		;field 8
		inv_lat = 0.
		READU, unit, inv_lat
		PRINTF, output_unit, "Invariant lat. = ", inv_lat
		;field 9
		mc_ilwain = 0.
		READU, unit, mc_ilwain
		PRINTF, output_unit, "Mc Ilwain param. = ", mc_ilwain
		;field 10
		geocentric_lat_conj = 0.
		READU, unit, geocentric_lat_conj
		PRINTF, output_unit, "Geocentric Lat. of the conjugate point = ", geocentric_lat_conj
		;field 11
		geocentric_long_conj = 0.
		READU, unit, geocentric_long_conj
		PRINTF, output_unit, "Geocentric Long. of the conjugate point = ", geocentric_long_conj
		;field 12
		geocentric_lat_conj_N_100 = 0.
		READU, unit, geocentric_lat_conj_N_100
		PRINTF, output_unit, "Geocentric Lat. of North conjugate point at 100km = ", $
			geocentric_lat_conj_N_100
		;field 13
		geocentric_long_conj_N_100 = 0.
		READU, unit, geocentric_long_conj_N_100
		PRINTF, output_unit, "Geocentric Long. of North conjugate point at 100km = ", $
			geocentric_long_conj_N_100
		;field 14
		geocentric_lat_conj_S_100 = 0.
		READU, unit, geocentric_lat_conj_S_100
		PRINTF, output_unit, "Geocentric Lat. of South conjugate point at 100km = ", $
			geocentric_lat_conj_S_100
		;field 15
		geocentric_long_conj_S_100 = 0.
		READU, unit, geocentric_long_conj_S_100
		PRINTF, output_unit, "Geocentric Long. of South conjugate point at 100km = ", $
			geocentric_long_conj_S_100
		;field 16
		B0 = FLTARR(3)
		READU, unit, B0
		PRINTF, output_unit, "B0 = ", B0
		;field 17
		Proton_geofrequency = 0.
		READU, unit, Proton_geofrequency
		PRINTF, output_unit, "Proton geofrequency = ", Proton_geofrequency
		;field 18
		Solar_pos = FLTARR(3)
		READU, unit, Solar_pos
		PRINTF, output_unit, "Solar position (Geog) = ", Solar_pos
		;field 19
		softv = 0UB
		READU, unit, softv
		;field 21
		softsv = 0UB
		READU, unit, softsv
		PRINTF, output_unit, "Soft version = ", STRING(FORMAT='(I1.1, ".", I1.1)', $
			softv, softsv)
	ENDIF ELSE BEGIN
		trash = BYTARR(90)
		READU, unit, trash
	ENDELSE

END
; End of read_block2
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block3
PRO read_block3, unit, output_unit, dh

	IF (dh EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block3"
		PRINTF, output_unit, "======"
		;field 1
		Msatgeog = FLTARR(3,3)
		READU, unit, Msatgeog
		PRINTF, output_unit, "Msatgeog = "
		PRINTF, output_unit, Msatgeog
		PRINTF, output_unit
		;field 2
		Mgeoggeom = FLTARR(3,3)
		READU, unit, Mgeoggeom
		PRINTF, output_unit, "Mgeoggeom = "
		PRINTF, output_unit, Mgeoggeom
		PRINTF, output_unit
		;Dummy = FLTARR(3,3)
		;READU, unit, Dummy
		;READU, unit, Dummy


		;field 5
		Quality = 0
		READU, unit, Quality
		PRINTF, output_unit, "Attitude quality = ", Quality
		;field 6
		softv = 0UB
		READU, unit, softv
		;field 7
		softsv = 0UB
		READU, unit, softsv
		PRINTF, output_unit, "Soft version = ", STRING(FORMAT='(I1.1, ".", I1.1)', $
			softv, softsv)
	ENDIF ELSE BEGIN
		trash = BYTARR(76)
		READU, unit, trash
	ENDELSE

END
; End of read_block3
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block4_ubf_wave
PRO read_block4_ubf_wave, unit, output_unit, dh, d

	IF (dh EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4"
		PRINTF, output_unit, "======"
		;field 1
		Data_type = BYTARR(21)
		READU, unit, Data_type
		Data_type = STRING(Data_type)
		PRINTF, output_unit, "Data type = ", Data_type
		;field 2
		Hks = BYTARR(32)
		READU, unit, Hks
		PRINTF, output_unit, Hks, FORMAT= '("HKs = ", 32Z2.2)'
		;field 3
		Coord_sys = BYTARR(9)
		READU, unit, Coord_sys
		Coord_sys = STRING(Coord_sys)
		PRINTF, output_unit, "Coordinate system = ", Coord_sys
		;field 4
		Msensat = FLTARR(3,3)
		READU, unit, Msensat
		PRINTF, output_unit, "Msensat = "
		PRINTF, output_unit, Msensat
		PRINTF, output_unit
		;field 5
		Data_unit = BYTARR(16)
		READU, unit, Data_unit
		Data_unit = STRING(Data_unit)
		PRINTF, output_unit, "Data unit = ", Data_unit
		;field 6
		Fe = 0. ; FLTARR(1)
		READU, unit, Fe
		PRINTF, output_unit, "Sampling frequency = ", Fe
		;field 7
		nbdata = 0 ; INTARR(1)
		READU, unit, nbdata
		PRINTF, output_unit, "Sample data per component = ", nbdata
		;field 8
		duration = 0. ; FLTARR(1)
		READU, unit, duration
		PRINTF, output_unit, "Time duration of one data array = ", duration
	ENDIF ELSE BEGIN
		trash = BYTARR(124)
		READU, unit, trash
	ENDELSE
	;field 9
	name1 = BYTARR(3)
	READU, unit, name1
	name1 = STRING(name1)
	;field 10
	data1 = FLTARR(256)
	READU, unit, data1
	;field 11
	name2 = BYTARR(3)
	READU, unit, name2
	name2 = STRING(name2)
	;field 12
	data2 = FLTARR(256)
	READU, unit, data2
	;field 13
	name3 = BYTARR(3)
	READU, unit, name3
	name3 = STRING(name3)
	;field 14
	data3 = FLTARR(256)
	READU, unit, data3
	;field 15
	name4 = BYTARR(3)
	READU, unit, name4
	name4 = STRING(name4)
	;field 16
	data4 = FLTARR(256)
	READU, unit, data4
	;field 17
	name5 = BYTARR(3)
	READU, unit, name5
	name5 = STRING(name5)
	;field 18
	data5 = FLTARR(256)
	READU, unit, data5
	;field 19
	name6 = BYTARR(3)
	READU, unit, name6
	name6 = STRING(name6)
	;field 20
	data6 = FLTARR(256)
	READU, unit, data6
	;field 21
	name7 = BYTARR(3)
	READU, unit, name7
	name7 = STRING(name7)
	;field 22
	data7 = FLTARR(256)
	READU, unit, data7
	IF (d EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4 data"
		PRINTF, output_unit, "======"
		PRINTF, output_unit, "Name 1 = ", name1
		PRINTF, output_unit, "data 1 = ", data1
		PRINTF, output_unit, "Name 2 = ", name2
		PRINTF, output_unit, "data 2 = ", data2
		PRINTF, output_unit, "Name 3 = ", name3
		PRINTF, output_unit, "data 3 = ", data3
		PRINTF, output_unit, "Name 4 = ", name4
		PRINTF, output_unit, "data 4 = ", data4
		PRINTF, output_unit, "Name 5 = ", name5
		PRINTF, output_unit, "data 5 = ", data5
		PRINTF, output_unit, "Name 6 = ", name6
		PRINTF, output_unit, "data 6 = ", data6
		PRINTF, output_unit, "Name 7 = ", name7
		PRINTF, output_unit, "data 7 = ", data7
	ENDIF $
	ELSE BEGIN
		IF (dh EQ 1) THEN BEGIN
			PRINTF, output_unit, "Name 1 = ", name1
			PRINTF, output_unit, "Name 2 = ", name2
			PRINTF, output_unit, "Name 3 = ", name3
			PRINTF, output_unit, "Name 4 = ", name4
			PRINTF, output_unit, "Name 5 = ", name5
			PRINTF, output_unit, "Name 6 = ", name6
			PRINTF, output_unit, "Name 7 = ", name7
		ENDIF
	ENDELSE

END
; End of read_block4_ubf_wave
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block4_ebf_wave
PRO read_block4_ebf_wave, unit, output_unit,dh, d

	IF (dh EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4"
		PRINTF, output_unit, "======"
		;field 1
		Data_type = BYTARR(21)
		READU, unit, Data_type
		Data_type = STRING(Data_type)
		PRINTF, output_unit, "Data type = ", Data_type
		;field 2
		Hks = BYTARR(32)
		READU, unit, Hks
		PRINTF, output_unit, Hks, FORMAT= '("HKs = ", 32Z2.2)'
		;field 3
		Coord_sys = BYTARR(9)
		READU, unit, Coord_sys
		Coord_sys = STRING(Coord_sys)
		PRINTF, output_unit, "Coordinate system = ", Coord_sys
		;field 4
		Msensat = FLTARR(3,3)
		READU, unit, Msensat
		PRINTF, output_unit, "Msensat = "
		PRINTF, output_unit, Msensat
		PRINTF, output_unit
		;field 5
		Data_unit = BYTARR(16)
		READU, unit, Data_unit
		Data_unit = STRING(Data_unit)
		PRINTF, output_unit, "Data unit = ", Data_unit
		;field 6
		Fe = 0.
		READU, unit, Fe
		PRINTF, output_unit, "Sampling frequency = ", Fe
		;field 7
		nbdata = 0
		READU, unit, nbdata
		PRINTF, output_unit, "Sample data per component = ", nbdata
		;field 8
		duration = 0.
		READU, unit, duration
		PRINTF, output_unit, "Time duration of one data array = ", duration
	ENDIF ELSE BEGIN
		trash = BYTARR(124)
		READU, unit, trash
	ENDELSE
	;field 9
	name1 = BYTARR(3)
	READU, unit, name1
	name1 = STRING(name1)
	;field 10
	data1 = FLTARR(4096)
	READU, unit, data1
	;field 11
	name2 = BYTARR(3)
	READU, unit, name2
	name2 = STRING(name2)
	;field 12
	data2 = FLTARR(4096)
	READU, unit, data2
	;field 13
	name3 = BYTARR(3)
	READU, unit, name3
	name3 = STRING(name3)
	;field 14
	data3 = FLTARR(4096)
	READU, unit, data3
	IF (d EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4 data"
		PRINTF, output_unit, "======"
		PRINTF, output_unit, "Name 1 = ", name1
		PRINTF, output_unit, "data 1 = ", data1
		PRINTF, output_unit, "Name 2 = ", name2
		PRINTF, output_unit, "data 2 = ", data2
		PRINTF, output_unit, "Name 3 = ", name3
		PRINTF, output_unit, "data 3 = ", data3
	ENDIF $
	ELSE BEGIN
		IF (dh EQ 1) THEN BEGIN
			PRINTF, output_unit, "Name 1 = ", name1
			PRINTF, output_unit, "Name 2 = ", name2
			PRINTF, output_unit, "Name 3 = ", name3
		ENDIF
	ENDELSE

END
; End of read_block4_ebf_wave
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block4_tbf_wave
PRO read_block4_tbf_wave, unit, output_unit, dh, d

	IF (dh EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4"
		PRINTF, output_unit, "======"
		;field 1
		Data_type = BYTARR(21)
		READU, unit, Data_type
		Data_type = STRING(Data_type)
		PRINTF, output_unit, "Data type = ", Data_type
		;field 2
		Hks = BYTARR(32)
		READU, unit, Hks
		PRINTF, output_unit, Hks, FORMAT= '("HKs = ", 32Z2.2)'
		;field 3
		Coord_sys = BYTARR(9)
		READU, unit, Coord_sys
		Coord_sys = STRING(Coord_sys)
		PRINTF, output_unit, "Coordinate system = ", Coord_sys
		;field 4
		Data_unit = BYTARR(16)
		READU, unit, Data_unit
		Data_unit = STRING(Data_unit)
		PRINTF, output_unit, "Data unit = ", Data_unit
		;field 5
		Fe = 0.
		READU, unit, Fe
		PRINTF, output_unit, "Sampling frequency = ", Fe
		;field 6
		nbdata = 0
		READU, unit, nbdata
		PRINTF, output_unit, "Sample data per component = ", nbdata
		;field 7
		duration = 0.
		READU, unit, duration
		PRINTF, output_unit, "Time duration of one data array = ", duration
	ENDIF ELSE BEGIN
		trash = BYTARR(88)
		READU, unit, trash
	ENDELSE
	;field 8
	name1 = BYTARR(3)
	READU, unit, name1
	name1 = STRING(name1)
	;field 9
	data1 = FLTARR(8192)
	READU, unit, data1
	IF (d EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4 data"
		PRINTF, output_unit, "======"
		PRINTF, output_unit, "Name 1 = ", name1
		PRINTF, output_unit, "data 1 = ", data1
	ENDIF $
	ELSE BEGIN
		IF (dh EQ 1) THEN BEGIN
			PRINTF, output_unit, "Name 1 = ", name1
		ENDIF
	ENDELSE

END
; End of read_block4_tbf_wave
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block4_spect
PRO read_block4_spect, unit, output_unit, dh, d

	;field 1
	Data_type = BYTARR(21)
	READU, unit, Data_type
	Data_type = STRING(Data_type)
	;field 2
	Hks = BYTARR(32)
	READU, unit, Hks
	;field 3
	Coord_sys = BYTARR(9)
	READU, unit, Coord_sys
	Coord_sys = STRING(Coord_sys)
	;field 4
	comp = BYTARR(3)
	READU, unit, comp
	comp = STRING(comp)
	;field 5
	Data_unit = BYTARR(16)
	READU, unit, Data_unit
	Data_unit = STRING(Data_unit)
	;field 6
	NbSpect = 0UB
	READU, unit, NbSpect
	;field 7
	NbComp = 0
	READU, unit, NbComp
	;field 8
	duration = 0.
	READU, unit, duration
	;field 9
	FreqRes = 0.
	READU, unit, FreqRes
	;field 10
	FreqRange = FLTARR(2)
	READU, unit, FreqRange
	;field 11
	Time = INTARR(7)
	READU, unit, Time
	IF (dh EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4"
		PRINTF, output_unit, "======"
		PRINTF, output_unit, "Data type = ", Data_type
		PRINTF, output_unit, Hks, FORMAT= '("HKs = ", 32Z2.2)'
		PRINTF, output_unit, "Coordinate system = ", Coord_sys
		PRINTF, output_unit, "Component name = ", comp
		PRINTF, output_unit, "Data unit = ", Data_unit
		PRINTF, output_unit, "Number of consecutice spectra = ", NbSpect
		PRINTF, output_unit, "Number of spectrum freauencies = ", NbComp
		PRINTF, output_unit, "Time duration of one data array = ", duration
		PRINTF, output_unit, "Frequency resolution = ", FreqRes
		PRINTF, output_unit, "Frequency range = ", FreqRange
		PRINTF, output_unit, "UT time = ", Time
	ENDIF
	;field 12
	Spectra = FLTARR(NbComp, NbSpect)
	READU, unit, Spectra
	IF (d EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4 data"
		PRINTF, output_unit, "======"
		FOR i = 0, NbSpect-1 DO BEGIN
			PRINTF, output_unit, "spectrum " + STRCOMPRESS(i, /REMOVE_ALL) + " = ", Spectra[*, i]
		ENDFOR
	ENDIF

END
; End of read_block4_spect
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block4_hf_wave
PRO read_block4_hf_wave, unit, output_unit, dh, d

	IF (dh EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4"
		PRINTF, output_unit, "======"
		;field 1
		Data_type = BYTARR(21)
		READU, unit, Data_type
		Data_type = STRING(Data_type)
		PRINTF, output_unit, "Data type = ", Data_type
		;field 2
		Hks = BYTARR(32)
		READU, unit, Hks
		PRINTF, output_unit, Hks, FORMAT= '("HKs = ", 32Z2.2)'
		;field 3
		Coord_sys = BYTARR(9)
		READU, unit, Coord_sys
		Coord_sys = STRING(Coord_sys)
		PRINTF, output_unit, "Coordinate system = ", Coord_sys
		;field 4
		Data_unit = BYTARR(16)
		READU, unit, Data_unit
		Data_unit = STRING(Data_unit)
		PRINTF, output_unit, "Data unit = ", Data_unit
		;field 5
		Fe = 0.
		READU, unit, Fe
		PRINTF, output_unit, "Sampling frequency = ", Fe
		;field 6
		nbdata = 0
		READU, unit, nbdata
		PRINTF, output_unit, "Sample data per component = ", nbdata
		;field 7
		duration = 0.
		READU, unit, duration
		PRINTF, output_unit, "Time duration of one data array = ", duration
	ENDIF ELSE BEGIN
		trash = BYTARR(88)
		READU, unit, trash
	ENDELSE
	;field 8
	name1 = BYTARR(3)
	READU, unit, name1
	name1 = STRING(name1)
	;field 9
	data1 = FLTARR(4096)
	READU, unit, data1
	IF (d EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4 data"
		PRINTF, output_unit, "======"
		PRINTF, output_unit, "Name 1 = ", name1
		PRINTF, output_unit, "data 1 = ", data1
	ENDIF $
	ELSE BEGIN
		IF (dh EQ 1) THEN BEGIN
			PRINTF, output_unit, "Name 1 = ", name1
		ENDIF
	ENDELSE
END
; End of read_block4_hf_wave
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block4_iap
PRO read_block4_iap, unit, output_unit, dh, d

	IF (dh EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4"
		PRINTF, output_unit, "======"
		;field 1
		Data_type = BYTARR(10)
		READU, unit, Data_type
		Data_type = STRING(Data_type)
		PRINTF, output_unit, "Data type = ", Data_type
		;field 2
		Hks = BYTARR(32)
		READU, unit, Hks
		PRINTF, output_unit, Hks, FORMAT= '("HKs = ", 32Z2.2)'
		;field 3
		Res = 0.
		READU, unit, Res
		PRINTF, output_unit, "Time resolution = ", Res
		;field 4
		Density_unit = BYTARR(6)
		READU, unit, Density_unit
		Density_unit = STRING(Density_unit)
		PRINTF, output_unit, "Density_unit = ", Density_unit
		;field 5
		Temperature_unit = BYTARR(6)
		READU, unit, Temperature_unit
		Temperature_unit = STRING(Temperature_unit)
		PRINTF, output_unit, "Temperature_unit = ", Temperature_unit
		;field 6
		Velocity_unit = BYTARR(6)
		READU, unit, Velocity_unit
		Velocity_unit = STRING(Velocity_unit)
		PRINTF, output_unit, "Velocity_unit = ", Velocity_unit
		;field 7
		Potential_unit = BYTARR(6)
		READU, unit, Potential_unit
		Potential_unit = STRING(Potential_unit)
		PRINTF, output_unit, "Potential_unit = ", Potential_unit
		;field 8
		Angle_unit = BYTARR(6)
		READU, unit, Angle_unit
		Angle_unit= STRING(Angle_unit)
		PRINTF, output_unit, "Angle_unit = ", Angle_unit
		;field 9
		H = 0.
		READU, unit, H
		PRINTF, output_unit, "H+ density = ", H
		;field 10
		He = 0.
		READU, unit, He
		PRINTF, output_unit, "He+ density = ", He
		;field 11
		O = 0.
		READU, unit, O
		PRINTF, output_unit, "O+ density = ", O
		;field 12
		IonTemp = 0.
		READU, unit, IonTemp
		PRINTF, output_unit, "Ions temperature = ", IonTemp
		;field 13
		IonVelocity = 0.
		READU, unit, IonVelocity
		PRINTF, output_unit, "Ions velocity along Oz axis = ", IonVelocity
		;field 14
		Angle_V_Oz = 0.
		READU, unit, Angle_V_Oz
		PRINTF, output_unit, "Angle (V, Oz) = ", Angle_V_Oz
		;field 15
		Angle_Vxy_Ox = 0.
		READU, unit, Angle_Vxy_Ox
		PRINTF, output_unit, "Angle (Vxy, Ox) = ", Angle_Vxy_Ox
		;field 16
		Vs = 0.
		READU, unit, Vs
		PRINTF, output_unit, "Satellite Potential = ", Vs
	ENDIF ELSE BEGIN
		trash = BYTARR(108)
		READU, unit, trash
	ENDELSE

END
; End of read_block4_iap
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block4_isl
PRO read_block4_isl, unit, output_unit, dh, d

  IF (dh EQ 1) THEN BEGIN
    PRINTF, output_unit, ""
    PRINTF, output_unit, "Block4"
    PRINTF, output_unit, "======"
                                ;field 1
    Data_type = BYTARR(10)
    READU, unit, Data_type
    Data_type = STRING(Data_type)
    PRINTF, output_unit, "Data type = ", Data_type
                                ;field 2
    Hks = BYTARR(32)
    READU, unit, Hks
    PRINTF, output_unit, Hks, FORMAT= '("HKs = ", 32Z2.2)'
                                ;field 3
    Res = 0.
    READU, unit, Res
    PRINTF, output_unit, "Time resolution = ", Res
                                ;field 4
    Density_unit = BYTARR(5)
    READU, unit, Density_unit
    Density_unit = STRING(Density_unit)
    PRINTF, output_unit, "Density_unit = ", Density_unit
                                ;field 5
    Temperature_unit = BYTARR(5)
    READU, unit, Temperature_unit
    Temperature_unit = STRING(Temperature_unit)
    PRINTF, output_unit, "Temperature_unit = ", Temperature_unit
                                ;field 6
    Potential_unit = BYTARR(5)
    READU, unit, Potential_unit
    Potential_unit = STRING(Potential_unit)
    PRINTF, output_unit, "Potential_unit = ", Potential_unit
                                ;field 7
  ENDIF ELSE BEGIN
    trash = BYTARR(61)
    READU, unit, trash
  ENDELSE
  IF (d EQ 1) THEN BEGIN
    Ed = 0.
    READU, unit, Ed
    PRINTF, output_unit, "Electron density = ", Ed
                                ;field 8
    Ni = 0.
    READU, unit, Ni
    PRINTF, output_unit, "Ion density = ", Ni
                                ;field 9
    Et = 0.
    READU, unit, Et
    PRINTF, output_unit, "Electron temperature = ", Et
                                ;field 10
    Vp = 0.
    READU, unit, Vp
    PRINTF, output_unit, "Plasma potential = ", Vp
                                ;field 11
    Vi0 = 0.
    READU, unit, Vi0
    PRINTF, output_unit, "Floating potential = ", Vi0
                                ;field 12
    Vs = 0.
    READU, unit, Vs
    PRINTF, output_unit, "Satellite potential = ", Vs
  ENDIF ELSE BEGIN
    trash = BYTARR(24)
    READU, unit, trash
  ENDELSE
  
END
; End of read_block4_isl
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block4_idp_burst
PRO read_block4_idp_burst, unit, output_unit, dh, d

	IF (dh EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4"
		PRINTF, output_unit, "======"
		;field 1
		Data_type = BYTARR(10)
		READU, unit, Data_type
		Data_type = STRING(Data_type)
		PRINTF, output_unit, "Data type = ", Data_type
		;field 2
		Hks = BYTARR(32)
		READU, unit, Hks
		PRINTF, output_unit, Hks, FORMAT= '("HKs = ", 32Z2.2)'
		;field 3
		Res = 0.
		READU, unit, Res
		PRINTF, output_unit, "Time resolution = ", Res
		;field 4
		Polarisation = 0.
		READU, unit, Polarisation
		PRINTF, output_unit, "Polarisation voltage (V) = ", Polarisation
		;field 5
		Discrimination = 0.
		READU, unit, Discrimination
		PRINTF, output_unit, "Discrimination level (keV) = ", Discrimination
		;field 6
		Spect_unit = BYTARR(20)
		READU, unit, Spect_unit
		Spect_unit = STRING(Spect_unit)
		PRINTF, output_unit, "Spectrum data unit = ", Spect_unit
		;field 7
		Pitch_angle_unit = BYTARR(6)
		READU, unit, Pitch_angle_unit
		Pitch_angle_unit = STRING(Pitch_angle_unit)
		PRINTF, output_unit, "Pitch angle unit = ", Pitch_angle_unit
	ENDIF ELSE BEGIN
		trash = BYTARR(80)
		READU, unit, trash
	ENDELSE
	;field 8
	data1 = FLTARR(256)
	READU, unit, data1
	;field 9
	data2 = FLTARR(256)
	READU, unit, data2
	;field 10
	data3 = FLTARR(256)
	READU, unit, data3
	;field 11
	data4 = FLTARR(256)
	READU, unit, data4
	;field 12
	Energy_table = FLTARR(256)
	READU, unit, Energy_table
	;field 13
	Pitch_angle = 0.
	READU, unit, Pitch_angle
	IF (d EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4 data"
		PRINTF, output_unit, "======"
		PRINTF, output_unit, "data 1 = ", data1
		PRINTF, output_unit, "data 2 = ", data2
		PRINTF, output_unit, "data 3 = ", data3
		PRINTF, output_unit, "data 4 = ", data4
		PRINTF, output_unit, "Energy_table = ", Energy_table
		PRINTF, output_unit, "Pitch angle = ", Pitch_angle
	ENDIF $
	ELSE BEGIN
		IF (dh EQ 1) THEN BEGIN
			PRINTF, output_unit, "Pitch angle = ", Pitch_angle
		ENDIF
	ENDELSE

END
; End of read_block4_idp_burst
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block4_idp_survey
PRO read_block4_idp_survey, unit, output_unit, dh, d

  IF (dh EQ 1) THEN BEGIN
    PRINTF, output_unit, ""
    PRINTF, output_unit, "Block4"
    PRINTF, output_unit, "======"
                                ;field 1
    Data_type = BYTARR(10)
    READU, unit, Data_type
    Data_type = STRING(Data_type)
    PRINTF, output_unit, "Data type = ", Data_type
                                ;field 2
    Hks = BYTARR(32)
    READU, unit, Hks
    PRINTF, output_unit, Hks, FORMAT= '("HKs = ", 32Z2.2)'
                                ;field 3
    ResS = 0.
    READU, unit, ResS
    PRINTF, output_unit, "Spectrum time resolution (s) = ", ResS
                                ;field 4
    ResC = 0.
    READU, unit, ResC
    PRINTF, output_unit, "Counters time resolution (s) = ", ResC
                                ;field 5
    Polarisation = 0.
    READU, unit, Polarisation
    PRINTF, output_unit, "Polarisation voltage (V) = ", Polarisation
                                ;field 6
    Discrimination = 0.
    READU, unit, Discrimination
    PRINTF, output_unit, "Discrimination level (keV) = ", Discrimination
                                ;field 7
    Threshold_low_1 = 0.
    READU, unit, Threshold_low_1
    PRINTF, output_unit, "Threshold low interval 1 (keV) = ", Threshold_low_1
                                ;field 8
    Threshold_low_2 = 0.
    READU, unit, Threshold_low_2
    PRINTF, output_unit, "Threshold low interval 2 (keV) = ", Threshold_low_2
                                ;field 9
    Threshold_low_3 = 0.
    READU, unit, Threshold_low_3
    PRINTF, output_unit, "Threshold low interval 3 (keV) = ", Threshold_low_3
                                ;field 10
    Threshold_high_3 = 0.
    READU, unit, Threshold_high_3
    PRINTF, output_unit, "Threshold high interval 3 (keV) = ", Threshold_high_3
                                ;field 11
    Spect_unit = BYTARR(20)
    READU, unit, Spect_unit
    Spect_unit = STRING(Spect_unit)
    PRINTF, output_unit, "Spectrum data unit = ", Spect_unit
                                ;field 12
    Pitch_angle_unit = BYTARR(6)
    READU, unit, Pitch_angle_unit
    Pitch_angle_unit = STRING(Pitch_angle_unit)
    PRINTF, output_unit, "Pitch angle unit = ", Pitch_angle_unit
  ENDIF ELSE BEGIN
    trash = BYTARR(100)
    READU, unit, trash
  ENDELSE
  
  IF (d EQ 1) THEN BEGIN
    PRINTF, output_unit, ""
    PRINTF, output_unit, "Block4 data"
    PRINTF, output_unit, "======"
    C1 = 0UL
    C2 = 0UL
    C3 = 0UL
    data = FLTARR(128)
    t = 0.
    FOR j= 0, 6 DO BEGIN
      ts = t
      PRINTF, output_unit, "      delta_t         Band1       Band2     Band3"
      PRINTF, output_unit, "      -------------------------------------------"
      FOR i = 0, 3 DO BEGIN     ;fields 13, 15, 17, 19, 21, 23, 25
        READU, unit, C1
        READU, unit, C2
        READU, unit, C3
        PRINTF, output_unit, t, C1, C2, C3
        t = t + 1.
      ENDFOR
      READU, unit, data         ;fields 14, 16, 18, 20, 22, 24, 26
      PRINTF, output_unit, "Spectrum at ", ts
      PRINTF, output_unit, data
    ENDFOR
                                ;field 27
    Energy_table = FLTARR(128)
    READU, unit, Energy_table
    PRINTF, output_unit, "Energy_table = ", Energy_table
                                ;field 28
    Pitch_angle = 0.
    READU, unit, Pitch_angle
    PRINTF, output_unit, "Pitch angle = ", Pitch_angle
  ENDIF $
  ELSE BEGIN
    IF (dh EQ 1) THEN BEGIN
      trash = BYTARR(4432)
      READU, unit, trash
      Pitch_angle = 0.
      READU, unit, Pitch_angle
      PRINTF, output_unit, "Pitch angle = ", Pitch_angle
    ENDIF $
    ELSE trash = BYTARR(4436)
    READU, unit, trash
  ENDELSE
  
END
; End of read_block4_idp_survey
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_block4_rnf
PRO read_block4_rnf, unit, output_unit, dh, d

	;field 1
	Data_type = BYTARR(21)
	READU, unit, Data_type
	Data_type = STRING(Data_type)
	;field 2
	Hks = BYTARR(32)
	READU, unit, Hks
	;field 3
	Sub = 0UB
	READU, unit, Sub
	;field 4
	Title = BYTARR(20)
	READU, unit, Title
	Title = STRING(Title)
	;field 5
	Comp = BYTARR(3)
	READU, unit, Comp
	Comp = STRING(Comp)
	;field 6
	Res = 0.
	READU, unit, Res
	;field 7
	Nbcl = 0UB
	READU, unit, Nbcl
	;field 8
	Nbs = 0UB
	READU, unit, Nbs
	;field 9
	Nbc = 0UB
	READU, unit, Nbc
	;field 10
	class = BYTARR(10)
	READU, unit, class
	class = STRING(class)
	;field 11
	Dmin = FLTARR(20)
	READU, unit, Dmin
	;field 12
	Dmax = FLTARR(20)
	READU, unit, Dmax
	IF (dh EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4"
		PRINTF, output_unit, "======"
		PRINTF, output_unit, "Data type = ", Data_type
		PRINTF, output_unit, Hks, FORMAT= '("HKs = ", 32Z2.2)'
		PRINTF, output_unit, "Data sub-type= ", Sub
		PRINTF, output_unit, "Study title = ", Title
		PRINTF, output_unit, "Component name = ", Comp
		PRINTF, output_unit, "Time resolution = ", Res
		PRINTF, output_unit, "Class number = ", Nbcl
		PRINTF, output_unit, 'Number spectra (or nb points when curves) = ', Nbs
		PRINTF, output_unit, '0 (or Curve number) = ', Nbc
		PRINTF, output_unit, "Unit name = ", class
		PRINTF, output_unit, "Range min for Di = ", Dmin
		PRINTF, output_unit, "Range max for Di = ", Dmax
	ENDIF
	;field 13
	Validity = BYTARR(128)
	READU, unit, Validity
	;field 14
	iCl = BYTARR(Nbcl, Nbs)
	READU, unit, iCl
	n = 2560L - LONG(Nbs) * LONG(Nbcl)
	Dummy = BYTARR(n)
	READU, unit, Dummy
	;field 15
	iCl2 = BYTARR(Nbcl, Nbs)
	READU, unit, iCl2
	IF (n GT 0) THEN READU, unit, Dummy
	IF (d EQ 1) THEN BEGIN
		PRINTF, output_unit, ""
		PRINTF, output_unit, "Block4 data"
		PRINTF, output_unit, "======"
		PRINTF, output_unit, "Validity Vector = "
		PRINTF, output_unit, Validity
		PRINTF, output_unit, "Spectrogram intensity = "
		PRINTF, output_unit, iCl
		PRINTF, output_unit, "Spectrogram uncertainty = "
		PRINTF, output_unit, iCl2
	ENDIF
END
; End of read_block4_rnf
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of read_n1
; Input(s) : apid : apid of the input file
; 			 unit : unit of the input file
; 			 output_unit : unit of the file where the result is printed
; 			 params : parameters corresponding to the blocks to print
; Output(s) :
PRO read_n1, apid, unit, output_unit, params

	POINT_LUN, unit, 0

	WHILE NOT EOF(unit) DO BEGIN
		read_block1, unit, t0, output_unit, params[0]
		read_block2, unit, output_unit, params[1]
		read_block3, unit, output_unit, params[2]
		CASE apid OF
			'1129' : read_block4_ubf_wave, unit, output_unit, params[3], params[4]
			'1130' : read_block4_ebf_wave, unit, output_unit, params[3], params[4]
			'1131' : read_block4_tbf_wave, unit, output_unit, params[3], params[4]
			'1132' : read_block4_spect, unit, output_unit, params[3], params[4]
			'1133' : read_block4_hf_wave, unit, output_unit, params[3], params[4]
			'1134' : read_block4_spect, unit, output_unit, params[3], params[4]
			'1135' : read_block4_ebf_wave, unit, output_unit, params[3], params[4]
			'1136' : read_block4_tbf_wave, unit, output_unit, params[3], params[4]
			'1137' : read_block4_spect, unit, output_unit, params[3], params[4]
			'1138' : read_block4_rnf, unit, output_unit, params[3], params[4]
			'1139' : read_block4_iap, unit, output_unit, params[3], params[4]
			'1140' : read_block4_iap, unit, output_unit, params[3], params[4]
                        '1141' : read_block4_idp_burst, unit, output_unit, params[3], params[4]
                        '1142' : read_block4_idp_survey, unit, output_unit, params[3], params[4]
                        '1143' : read_block4_isl, unit, output_unit, params[3], params[4]
                        '1144' : read_block4_isl, unit, output_unit, params[3], params[4]
                      ENDCASE
	ENDWHILE

END
; End of read_n1
;------------------------------------------------------------------------------


;------------------------------------------------------------------------------
; Begin of print_n1
; Input(s) : st : structure defined in rd_dmt_n1
; 			 output_unit : unit of the temporary file where the result is printed
; Output(s) :
PRO print_n1, st, output_unit

	params = st.val_param
	packet_nb = st.val_packet
	apid = st.apid

	res = FSTAT(st.file_input_unit)
	cur_fp = 0L ; cur_fp like 'CURrent position of the File Pointer'
	POINT_LUN, st.file_input_unit, cur_fp
	i = 0L

	WHILE (EOF(st.file_input_unit) EQ 0) AND (i LT packet_nb) DO BEGIN
		; Block 1, 2 and 3
		old_fp = cur_fp
		cur_fp += 204L
		POINT_LUN, st.file_input_unit, cur_fp
		; Block 4
		CASE apid OF
			'1129' : cur_fp += 7313L
			'1130' : cur_fp += 49285L
			'1131' : cur_fp += 32859L
			'1132' : BEGIN
				cur_fp += 81L
				POINT_LUN, st.file_input_unit, cur_fp
				nb = 0B
				READU, st.file_input_unit, nb
				nbf = 0
				READU, st.file_input_unit, nbf
				; 3 : nb and nbf, 30 : fields 8 to 11, nbf*4*nb : power spectrum data
				cur_fp += (3L + 30L + LONG(Nbf)*4L*LONG(nb))
			END
			'1133' : cur_fp += 16475L
			'1134' : BEGIN
				cur_fp += 81L
				POINT_LUN, st.file_input_unit, cur_fp
				nb = 0B
				READU, st.file_input_unit, nb
				nbf = 0
				READU, st.file_input_unit, nbf
				; 3 : nb and nbf, 30 : fields 8 to 11, nbf*4*nb : power spectrum data
				cur_fp += (3L + 30L + LONG(Nbf)*4L*LONG(nb))
			END
			'1135' : cur_fp += 49285L
			'1136' : cur_fp += 32859L
			'1137' : BEGIN
				cur_fp += 81L
				POINT_LUN, st.file_input_unit, cur_fp
				nb = 0B
				READU, st.file_input_unit, nb
				nbf = 0
				READU, st.file_input_unit, nbf
				; 3 : nb and nbf, 30 : fields 8 to 11, nbf*4*nb : power spectrum data
				cur_fp += (3L + 30L + LONG(Nbf)*4L*LONG(nb))
			END
			'1138' : cur_fp += 5502L
			'1139' : cur_fp += 108L
			'1140' : cur_fp += 108L
			'1141' : cur_fp += 5204L
			'1142' : cur_fp += 4536L
			'1143' : cur_fp += 85L
			'1144' : cur_fp += 85L
		ENDCASE
		POINT_LUN, st.file_input_unit, cur_fp
		i += 1L
	ENDWHILE
	IF (EOF(st.file_input_unit) EQ 1) THEN BEGIN
		st.val_packet = i-1L
		POINT_LUN, st.file_input_unit, old_fp
	ENDIF ELSE st.val_packet = i
	read_block1, st.file_input_unit, t0, output_unit, params[0]
	read_block2, st.file_input_unit, output_unit, params[1]
	read_block3, st.file_input_unit, output_unit, params[2]

        CASE apid OF
		'1129' : read_block4_ubf_wave, st.file_input_unit, output_unit, params[3], params[4]
		'1130' : read_block4_ebf_wave, st.file_input_unit, output_unit, params[3], params[4]
		'1131' : read_block4_tbf_wave, st.file_input_unit, output_unit, params[3], params[4]
		'1132' : read_block4_spect, st.file_input_unit, output_unit, params[3], params[4]
		'1133' : read_block4_hf_wave, st.file_input_unit, output_unit, params[3], params[4]
		'1134' : read_block4_spect, st.file_input_unit, output_unit, params[3], params[4]
		'1135' : read_block4_ebf_wave, st.file_input_unit, output_unit, params[3], params[4]
		'1136' : read_block4_tbf_wave, st.file_input_unit, output_unit, params[3], params[4]
		'1137' : read_block4_spect, st.file_input_unit, output_unit, params[3], params[4]
		'1138' : read_block4_rnf, st.file_input_unit, output_unit, params[3], params[4]
		'1139' : read_block4_iap, st.file_input_unit, output_unit, params[3], params[4]
		'1140' : read_block4_iap, st.file_input_unit, output_unit, params[3], params[4]
		'1141' : read_block4_idp_burst, st.file_input_unit, output_unit, params[3], params[4]
		'1142' : read_block4_idp_survey, st.file_input_unit, output_unit, params[3], params[4]
		'1143' : read_block4_isl, st.file_input_unit, output_unit, params[3], params[4]
		'1144' : read_block4_isl, st.file_input_unit, output_unit, params[3], params[4]
	ENDCASE

END
; End of print_n1
;------------------------------------------------------------------------------


PRO rd_dmt_n1_sub
END
