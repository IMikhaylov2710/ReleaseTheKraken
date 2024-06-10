import os


#apt-get install build-essential autoconf automake libtool

curDir = os.path.dirname(os.path.realpath(__file__))
#kraken2
os.system('cd %s && git clone https://github.com/DerrickWood/kraken2/' % curDir)
krakenDir = os.path.join(curDir, 'kraken2')
os.system('cd %s && ./install_kraken2.sh .' % krakenDir)
print('cd %s && kraken2-build --db silva --special silva' % krakenDir)
os.system('cd %s && kraken2-build --db silva --special silva' % krakenDir)