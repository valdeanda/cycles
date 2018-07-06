#
#
#
# mie ene 22 09:05:56 CST 2014
#

import os
import re
import urllib2
import numpy as np
import matplotlib.pyplot as plt
import pandas as  pd


def actualize_org_list():
    """ Actualize from KEGG the file that contains the list of organisms
    containded in the database
    """
    if not os.path.exists('./organism.txt'):
        question = 'y'
    else:
        question = raw_input("Overwrite the previous organism.txt file (y/n): ")
    if question == 'y':
        print "Downloading and storing species list..."
        with open('organism.txt', 'w') as outf:
            handler = urllib2.urlopen("http://rest.kegg.jp/list/organism")
            outf.write(handler.read())
        print "--- Done!!!\n "
    elif question == 'n':
        print "The organism.txt was not actualized."
    else:
        print "Bad answer -> {}.".format(question)

def actualize_ec_list():
    """Actualize from KEGG the file that contains the list of EC numbers
    containded in the database
    """
    if not os.path.exists('./ec.txt'):
        question = 'y'
    else:
        question = raw_input("Overwrite the previous ec.txt file (y/n): ")
    if question == 'y':
        print "Downloading and storing EC number list..."
        with open('ec.txt', 'w') as outf:
            handler = urllib2.urlopen("http://rest.kegg.jp/list/ec")
            outf.write(handler.read())
        print "--- Done!!!\n "
    elif question == 'n':
        print "The ec.txt was not actualized."
    else:
        print "Bad answer -> {}.".format(question)
        
def parse__organism():
    """Returns a dictionary describing the data contained in
    organism.txt file. Keys = columns ; value = list of values
    """
    if not os.path.exists("organism.txt"):
        print "The organism.txt does not exists"
    else:
        data = {'id':[], 'Domain':[], 'Name':[], 'Division':[], "Subdivision":[]}
        with open('organism.txt') as inf:
            for line in inf:
                line = line.strip()
                elements = line.split('\t')
                data['id'].append(elements[1])
                data['Name'].append(elements[2])
                tax = elements[3].split(';')
                if tax[0] == 'Eukaryotes':
                    data['Domain'].append(tax[0])
                    data['Division'].append(tax[1])
                    data['Subdivision'].append(tax[2])
                elif tax[0] == 'Prokaryotes':
                    data['Domain'].append(tax[1])
                    data['Division'].append(tax[2])
                    try:
                        data['Subdivision'].append(tax[3])
                    except:
                        data['Subdivision'].append("Unclassified")
                else:
                    print "Oops!! some thing is werid"
    data = pd.DataFrame(data)
    return data
        
def ec_profile_by_ec(ec_list, actualize=False, save_files=True):
    """Returns a pandas DataFrame that contains the profiles of
    precense absence of the EC numbers in ec_list.
    
    Arguments:
    - `ec_list`:
    - `asctualize`: Bool, default False. If True, actualize the organism.txt
    file. CAUTION!! if actualize and save filesthe data will be overwriten
    - `save_files`: Bool, default False. If True, save the files downloaded
    by kegg
    """
    # load organism data and actualize if requiered
    if actualize:
        print "Actualizing organism.list"
        actualize_org_list()
        organisms = parse__organism()
    else:
        assert os.path.exists("organism.txt"), "organism.txt file does not exist\n"\
            "Run again with  actualize=True"
        organisms = parse__organism()
    # create DataFrame
    frame = pd.DataFrame(np.zeros((len(organisms), len(ec_list))),
                                  columns=ec_list,  dtype=np.bool)
    data = pd.concat([organisms, frame], axis=1)
    data.set_index('id', inplace=True)
    # compile the RE to find sp code
    spre = re.compile(r"\t([a-z]+):")
    #
    for ec in ec_list:
        file_path = './Data/by_ec/ec_{}.txt'.format(ec)
        if os.path.exists(file_path):
            if actualize:
                print "If save_files == True "\
                      "the data files will be overwriten\n--"
                print "Retriving EC:{} data ...".format(ec)
                url = "http://rest.kegg.jp/link/genes/ec:{}".format(ec)
                handler = urllib2.urlopen(url)
                text = handler.read()
                if save_files:
                    with open(file_path, 'w') as outf:
                        outf.write(text)
                    print "---- {} saved!".format(file_path)
            else:
                print "Loading file {}...".format(file_path)
                with open(file_path) as inf:
                    text = inf.read()
        else:
            print "Retriving EC:{} data ...".format(ec)
            url = "http://rest.kegg.jp/link/genes/ec:{}".format(ec)
            handler = urllib2.urlopen(url)
            text = handler.read()
            if save_files:
                with open(file_path, 'w') as outf:
                    outf.write(text)
                print "---- {} saved!".format(file_path)
        sps = set(spre.findall(text))
        for sp in sps:
            data.ix[sp, ec] = True
    # rearrange indexes
    data.reset_index(inplace=True)
    data.set_index(['Domain', 'Division', 'Subdivision', 'Name','id'],
                   inplace=True)
    return data

def plot_heatmap(dataFrame, edgecolor='w', cmap=plt.cm.jet):
    """Plot the heatmap of a data frame
    
    Arguments:
    - `DataFrame`:  pandas DataFrame
    - `edgecolor`: If 'none', no edge color.
    """
    ax = plt.pcolor(dataFrame, cmap=cmap, edgecolor=edgecolor)
    xlen = len(dataFrame.ix[0])
    ylen = len(dataFrame)
    plt.xlim(0, xlen)
    plt.ylim(0, ylen)    
    plt.xticks(np.arange(xlen)+0.5, dataFrame.columns, rotation=90)
    plt.yticks(np.arange(ylen)+0.5, dataFrame.index)
    ax = plt.gca()
    ax.invert_yaxis()
    cbar = plt.colorbar(shrink=0.5, aspect=15)
    plt.draw()
    return ax, cbar

def main():
    """Main function
    """
    datafolder = './Data/'
    if not os.path.exists(datafolder):
        os.mkdir(datafolder)
    ecdata = './Data/by_ec/'
    if not os.path.exists(ecdata):
        os.mkdir(ecdata)

    orgdata = './Data/by_sp/'
    

    
    
if __name__ == '__main__':
    main()
    ec_list = ["1.1.1.1", "1.1.1.100", "1.1.1.28", "1.1.1.305", "2.1.2.2", "2.10.1.1" ]
    val = ("1.1.1.310", "1.1.1.313", "1.1.1.337","1.1.1.338","1.12.98.4","1.13.11.18","1.13.11.55","1.14.11.17","1.14.13.111","1.14.13.131", "1.14.14.5", "1.2.1.73", "1.2.1.81","1.4.99.2", "1.5.1.38","1.8.1.17", "1.8.1.2", "1.8.2.1", "1.8.2.2","1.8.2.3", "1.8.2.4", "1.8.3.4", "1.8.4.10", "1.8.4.8", "1.8.5.2","1.8.5.3", "1.8.5.4", "1.8.7.1", "1.8.99.1","1.8.99.2", "1.8.99.3", "2.1.1.107", "2.1.1.251", "2.1.1.269", "2.3.1.8","2.3.1.30", "2.3.1.46", "2.3.3.15", "2.5.1.6", "2.5.1.47", "2.5.1.48","2.6.1.55","2.6.1.77", "2.7.1.25", "2.7.7.4","2.7.7.5", "2.8.1.1", "2.8.1.5", "2.8.1.3", "3.12.1.1","3.6.2.1", "3.5.5.8", "4.1.1.79","4.4.1.24","4.4.1.25", "4.4.1.3")
    
    ec_re = re.compile(r"ec:(\d+.\d+.\d+.\d+)\t")
    ec_file = open('ec.txt')
    ec_text = ec_file.read()
    ec_file.close()
    ec_list = ec_re.findall(ec_text)
    
    def download_all_ec(ec_list):
        """Download all ec number to gene relation
        
        Arguments:
        - `ec_list`: list of ec numbers
        """
        copy_list = ec_list[:]
        for ec in ec_list:
            file_path = './Data/by_ec/ec_{}.txt'.format(ec)
            if os.path.exists(file_path):
                print ">{}: Alredy exists, skip!".format(ec)
                continue
            print "Retriving EC:{} data ...".format(ec)
            try:
                url = "http://rest.kegg.jp/link/genes/ec:{}".format(ec)
                handler = urllib2.urlopen(url)
                text = handler.read()
            except:
                raise
            with open(file_path, 'w') as outf:
                outf.write(text)
            print "... {} saved!".format(file_path)


