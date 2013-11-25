"""
    ISO Parsing
    
    Reads data in from the .txt files from Andy's ISO Data.  Either 
     - outputs them in more convenient formats for [R]
    - converts them to objects for the optimization problem
"""
import csv, pdb, sys, collections

def getFuelType( gen_name ):
    """Translates: AES LOND_GR_RIDGE_9_99Steam to ("AES LOND_GR_RIDGE_9_99", "Steam") """
    types = ("Steam", "CT", "Diesel", "Hydro", "Nuclear", "Fixed_Import")
    for type in types:
        if type in gen_name:
            fuel_type = type
            indx = gen_name.find(type)
            gen_name = gen_name[:indx]
            break
    else: #matches the for
        fuel_type = "Other"

    return (gen_name, fuel_type)    

def bidInfoToDict( size_path, price_path):
    """Converts the GAMS EnergyBidSize.txt file and pEngeryBidPrice.txt files to a dictionary
    {generator: { time: [(block_num, block_size, block_price)] } }
    """
    #this is done in a fairly inefficient way, but who cares??
    size_file = open(size_path, 'rU')
    #burn a line.
    size_file.readline()

    bid_dict = {}
    #First go thorugh and create the generator entry and a block entry
    for line in size_file:
        #end of file signified by this character
        if "/;" in line:
            break

        #the following is seemingly the easiest way to split the line, even though it seems assinine.
        indx = line.rfind(" ")
        gen_name, block_size = line[:indx], line[indx:]
        block_size = float( block_size.strip() )

        #split the long string and clean it up
        gen_name, block_num, time =gen_name.strip().split(".")
        gen_name, fuel_type = getFuelType(gen_name)
        time = time.strip()
        
        #add internal dictionary and check for duplicates
        if gen_name not in bid_dict:
            bid_dict[gen_name] = {}

        if time in bid_dict[gen_name]:
            #make sure that this block has not yet been added
            t = filter( lambda (num, size, price) : num==block_num, bid_dict[gen_name][time] )
            assert not t
        else:
            bid_dict[gen_name][time] = []

        #we use None to represent no price data for now
        bid_dict[gen_name][time].append ( (block_num, block_size, None)  )

    #go back through add price data for everyone. 
    price_file = open(price_path, 'rU')
    price_file.readline()  #burn a line

    for line in price_file:
        #end of file signified by this character
        if "/;" in line:
            break

        #the following is seemingly the easiest way to split the line, even though it seems assinine.
        indx = line.rfind(" ")
        gen_name, block_price = line[:indx], line[indx:]
        block_price = float( block_price.strip() )

        #split the long string and clean it up
        gen_name, block_num, time =gen_name.strip().split(".")
        gen_name, fuel_type = getFuelType(gen_name)
        time = time.strip()

        if gen_name not in bid_dict:
            raise ValueError("Generator %s not found in price-data" %gen_name)
        if time not in bid_dict[gen_name]:
            raise ValueError("Time %s not found for generator %s" % time, gen_name)

        for ix, (iBlock, iSize, iPrice) in enumerate(bid_dict[gen_name][time]):
            if iBlock==block_num:
                if iPrice is not None:
                    raise ValueError("Multiple prices")
                else:
                    bid_dict[gen_name][time][ix] = (iBlock, iSize, block_price)
                    break
        else:
            raise ValueError("Block %s not found for time %s generator %s" % (iBlock, time, gen_name) )

    #confirm that everyone actually has a price
    for gen, time_dicts in bid_dict.items():
        for bid_list in time_dicts.values():
            for (num, size, price) in bid_list:
                if price is None:
                    print "WTF!"
    
    return bid_dict

def readLoads( load_path ):
    """Reads the path pRandomAverage_v2.txt  """
    file = open(load_path, 'rU')
    #burn a line.
    file.readline()
    load_counter = collections.Counter()
    for line in file:
        #end of file signified by this character
        if "/;" in line:
            break
        entries = line.split()
        assert len(entries) == 2
        id = entries[0]
        indx = id.rfind(".") + 2
        load_counter[ int(id[indx:]) - 1 ] += float(entries[1]) 
    return load_counter

if __name__ == "__main__":
    pdb.set_trace()
    bid_dict = bidInfoToDict("../ISO-NE/pEnergyBidSize.txt", "../ISO-NE/pEnergyBidPrice.txt")

    print "Num Generators:\t", len(bid_dict)

    load_counter = readLoads("../ISO-NE/pRandomAverage_v2.txt")    
    print "Load Info:"
    for ix in xrange(1, 37):
        print ix, load_counter[ix]
    