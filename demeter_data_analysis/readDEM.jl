# julia script to read DEMETER_data
import Dates

# ----------------------------------------------------------------------------------------------
# function to open the files and read in packets into dict structure
function get_dataDict(datafile)
    println(datafile)
    f = open(datafile)
    lines = readlines(f)

    inds = [15:1:33;]
    data = Dict()

    for p in 1:1429:size(lines)[1] 
        n = convert(Int64, p ÷ 1429)
        packetn = "packet " * string(n)
        data[packetn] = Dict()

        # grab the header info that might be useful!
        for ind in inds
            currentline = lines[ind+p]
            keyarr = split(currentline,"=")
            # define the value
            if ind == 30 || ind == 32
                val = map(y -> parse(Float64,y), split(keyarr[2]))
            else 
                val = parse(Float64, keyarr[2])
            end
            keystr = keyarr[1] # define the string 
            push!(data[packetn], keystr => val) 
        end

        # add in the time, sampling freq, and unit
        time = map(y -> parse(Int64,y), split(split(lines[6+p],"e")[2]))
        timeobj = Dates.DateTime(time[1],time[2],time[3],time[4],time[5],time[6],time[7])
        data[packetn]["time"] = timeobj
        data[packetn]["fs"] = 40e3 # Hz 
        data[packetn]["unit"] = "uV/m"
        
        # add in the time domain data
        packetdata = lines[p+63:p+1428]
        packetdata[1] = packetdata[1][9:end]
        packetdatacleaned = collect(Iterators.flatten(map(y-> split(y), packetdata)))
        packetdatafinal = map(y-> parse(Float64, y), packetdatacleaned)   
        data[packetn]["Edata"] = packetdatafinal
    end
    return data
end
# ----------------------------------------------------------------------------------------------

# specify where the data is located
datadir = "/home/rileyannereid/workspace/canvas/DEMETER_data/data/inan_data/"

# get all the files in the directory
allfiles = String[]
for (root, dirs, files) in walkdir(datadir)
    append!(allfiles,joinpath.(root,files))
end

# check which have been parsed
datafiles = String[]
map(f -> if f[end-6:end-4] == "lec" push!(datafiles, f) end, allfiles)

@time get_dataDict(datafiles[1])