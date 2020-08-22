from Bio import SeqIO
import re

#Estrutura de Dados Trie - Busca a maior sequência conservada de aminoácios em cada site já identificado
class trieNo(object):

    def __init__(self, chtr):
        self.chr = chtr
        self.seq = []
        self.fimpalavra = False
        self.contaSeq = 1
        self.contaAmino = 1
        self.porcenCon = 60
        self.tipo = 3
        self.ranking = []

    #lê a sequencia de aminos no arquivo Fasta
    def leSequencias(self):
        listaSequencias = []
        for i in SeqIO.parse('seq-prot-ip3-R3.fas', 'fasta'):
            listaSequencias.append(list(i.seq))
        return listaSequencias


    #formata e adiciona os nomes das espécies à uma lista
    def leEspecies(self):
        listaEspecies = []
        for i in SeqIO.parse('seq-prot-ip3-R3.fas', 'fasta'):
            loc = re.search(r'type_3', i.id).end()
            sp = i.id[loc + 1:len(i.id)]
            sp = sp.split('_')
            sep = ' '
            item = sep.join(sp)
            listaEspecies.append(item)

        return listaEspecies


    def contaEspecies(self):
        listEsp = len(self.leEspecies())
        self.ranking = [0]*listEsp
        return listEsp



    # porcentagem a se considerar como site conservado - 60%
    def porcConservacao(self, numeroEspecies):
        return int(numeroEspecies * self.porcenCon/100)


    # elimina sequências com baixa conservação
    def locaisConserv(self, porConserv):
        listaSequencias = []
        for i in SeqIO.parse('seq-prot-ip3-R3.fas', 'fasta'):
            listaSequencias.append(i.seq)
        locaisConservados = []
        sitesConservados = []

        for j in range(len(listaSequencias[0])):
            linhas = []
            for k in range(len(listaSequencias)):
                if (listaSequencias[k][j] != '-'):
                    linhas.append(listaSequencias[k][j])

            if len(linhas) >= porConserv:
                sitesConservados.append(linhas.copy())
                locaisConservados.append(j)
            linhas.clear()
        return locaisConservados


    # agrupa os sites conservados consecutivos
    def agrupaSitesCons(self, listCon):
        listaTemp = []
        listaFim = []
        m = 0
        incr = 0
        tam = len(listCon)

        while incr < tam:
            chave = incr
            while listCon[incr] - listCon[chave] == m:
                listaTemp.append(listCon[incr])
                m = m + 1
                incr = incr + 1
                if incr == tam:
                    break
            listaFim.append(listaTemp.copy())
            listaTemp.clear()
            m = 0
        return listaFim


    # retorna uma lista com os sites conservados por espécie, uma linha por espécie
    def listAminosConservados(self, lisSitesCons, lisSeqAminos):
        listaTemp = []
        listaFim = []

        for i in range(len(lisSeqAminos)):
            for j in range(len(lisSeqAminos[i])):
                if j in lisSitesCons and lisSeqAminos[i][j] != '-':
                    listaTemp.append(lisSeqAminos[i][j])
            listaFim.append(listaTemp.copy())
            listaTemp.clear()
        return listaFim

    def maiorSequenciaConservada(self, lisCons, listSeq):
        maior = 0
        indice = 0
        listaFim = []
        for i in range(len(lisCons)):
            if len(lisCons[i]) > maior:
                maior = len(lisCons[i])
                indice = i

        inicio = lisCons[indice][0]
        final = lisCons[indice][maior-1]

        indNaoVazio = 0

        for j in range(len(listSeq)):
            if listSeq[j][inicio] != '-':
                indNaoVazio = j
                break
        listaFim.append(listSeq[indNaoVazio][inicio:final+1])
        sep = ''
        saida = sep.join(listaFim[0])


        return [maior,inicio, final, saida]

    def adicionaNo(self, palavra):
        node = self
        for caract in palavra:
            achou = False
            for item in node.seq:
                if item.chr == caract:
                    item.contaSeq = item.contaSeq+1
                    node = item
                    achou = True
                    break
            if not achou:
                novoNo = trieNo(caract)
                node.seq.append(novoNo)
                node = novoNo
        node.fimpalavra = True
        return [node.contaSeq, palavra]


    #pega os aminoacidos conservados por linha dentro da taxa de conservação estabelecida:
    def geraConserv(self, listaCon, listaSeq):
        listTemp = []
        listaAdd = []


        for item in listaCon:
            min = item[0]
            tam = len(item) -1
            max = item[tam]
            for p in range(len(listaSeq)):
                listTemp.append(listaSeq[p][min:max+1])

            listaAdd.append(listTemp.copy())
            listTemp.clear()


        #converte cada linha em uma palavra pra ser usada na Trie
        lTemp = []
        listaStr = []
        separator = ''
        for item in listaAdd:
            for j in item:
                strAm = separator.join(j)
                lTemp.append(strAm)
            listaStr.append(lTemp.copy())
            lTemp.clear()

        return listaStr



#gera as sequências mais encontradas em cada grupo de sites
    def geraPred(self, listStr):

        seqPred = ''
        listaSaidaPred = []

        for item in listStr:
            maior = 1
            for j in item:
                value = trieNo.adicionaNo(self, j)[0]
                nulo = j.find('-')
                if value > maior and nulo <= -1:
                    maior = value
                seqPred = trieNo.adicionaNo(self, j)[1]
            listaSaidaPred.append(seqPred)

        return listaSaidaPred


#gera ranking de espécies mais conservadas

    def insereAlte(self, indice):
        val = self.ranking[indice]
        val = val + 1
        self.ranking[indice] = val





    def geraRankOrd(self, listaEsp):

        vetEspOrd = [i for i in listaEsp]
        vetEspDes = [i for i in listaEsp]


        vetOrd = [int(i) for i in self.ranking]
        vetDes = [int(i) for i in self.ranking]

        for i in range(len(vetOrd)):

            min = i

            for j in range(i + 1, len(vetOrd)):
                if vetOrd[min] > vetOrd[j]:
                    min = j

            temp = vetOrd[i]
            vetOrd[i] = vetOrd[min]
            vetOrd[min] = temp

            temp2 = vetEspOrd[i]
            vetEspOrd[i] = vetEspOrd[min]
            vetEspOrd[min] = temp2




        return [vetDes, vetOrd, vetEspDes, vetEspOrd]






    def compara(raiz, listPred, listAminos, listEsp):

        node = trieNo(listPred)
        node.adicionaNo(listPred)

        listaAltFim = []
        listAminAlt = []
        listEspAlt = []
        listEspCons = []
        listNumAltAmin = []
        indEsp = 0

        alterSites = 0
        alterAmino = 0


        for j in range(len(listAminos)):
            values = listAminos[j]
            alterado = False
            result = []
            for letra in values:

                for item in node.seq:
                    if item.chr[0] == letra:
                        result.append(letra)
                    else:
                        result.append('[{}]'.format(letra))
                        alterado = True
                        alterAmino = alterAmino + 1
                    node = item


            if alterado:
                listAminAlt.append(result)
                alterSites = alterSites + 1
                listNumAltAmin.append(alterAmino)
                alterAmino = 0
                if result:
                    listaAltFim.append(result.copy())
                    listEspAlt.append(listEsp[j])
                    raiz.insereAlte(j)


            else:
                listEspCons.append(listEsp[j])


            node.adicionaNo(listPred)



            listAminAlt.clear()
            result.clear()


        node = None

        return [listaAltFim, listEspAlt, listNumAltAmin, alterSites, listEspCons]







root = trieNo('+')

especiesList = root.leEspecies()
qtdeEspecies = root.contaEspecies()
porceCon = root.porcConservacao(qtdeEspecies)
locaisConservados = root.locaisConserv(porceCon)
seqAminos = root.leSequencias()
indices = root.agrupaSitesCons(locaisConservados)
conservados = root.geraConserv(indices, seqAminos)
predominantes = root.geraPred(conservados)




# gera arquivo de saida

arquivo = "saidaIndels - Tipo"+str(root.tipo)+".txt"
with open(arquivo, "w") as arq:


    text09 = '\nSequências Conservadas, Inserções e Deleções\n'

    arq.write(text09)


    for j in range(len(indices)) :
        arq.write('\nSites:\n')
        arq.write(str(indices[j]))
        arq.write('\n')
        arq.write('\nSequência conservada mais comum:\n')
        arq.write(predominantes[j])
        arq.write('\n')





        arq.write('\nResumo do alinhamento na sequência:\n')
        listConser = root.compara(predominantes[j], conservados[j], especiesList)[4]
        arq.write('\n - Espécies Conservadas:\n')
        for item in range(len(listConser)):
            cons = str(listConser[item])
            arq.write(cons)
            arq.write('\n')

        valSitAlt = root.compara(predominantes[j], conservados[j], especiesList)[3]

        arq.write('\n - Espécies com INDELs:\n')
        arq.write('\nTotal:\n')
        arq.write(str(valSitAlt))
        arq.write('\n\n')




        list1 = root.compara(predominantes[j], conservados[j], especiesList)[0]
        list2 = root.compara(predominantes[j], conservados[j], especiesList)[1]
        list3 = root.compara(predominantes[j], conservados[j], especiesList)[2]


        for item in range(len(list1)):

            seq =  str(list1[item])
            esp = str(list2[item])
            alt = str(list3[item])


            arq.write('Espécie: {} -  Seq..: {} - INDELs: {}'.format(esp, seq, alt))

            arq.write('\n')

        arq.write('\n')


        arq.write('--------------------------------------------------------------------------------------------------------------------------')
        arq.write('\n\n')



    arq.close()













"""
Fluxo:

- Lê sequências 
- Lê espécies
- Conta quantidade de espécies
- Define porcentagem de conservação a considerar - porConserv
- Calcula a maior sequência conservada - maiorSequenciaConservada(lisCons)
- Elimine da lista os sites com baixa conservação - locaisConserv(porConserv)
- Sequencia os sites conservados - agrupaSitesCons(listCon):
- Gera lista com os aminos para cada grupo de sites conservados sequnciados anteriormente e Concatena os aminos, 
  convertendo em palavras, para serem usados na Trie - - geraConserv(listaCon, listaSeq):
  
- Insere a lista com os nomes na função adicionaNo para busca da maior sequência - geraPred(listaStr)
- Compara lista gerada de aminos nos sites conservados com a sequência predominate - compara(listaAminos, seqPred)
- Gerar arquivo de saida


"""
















