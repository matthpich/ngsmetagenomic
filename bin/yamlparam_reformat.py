import sys

with open(sys.argv[1], "r") as fi:
    output_l = []
    for line in fi:
        line = line.strip()
        if len(line) == 0: continue
        items = line.split(":", 1)
        items_ = []
        for item in items:
            item = item.strip()
            try:
                int(item)
                if item[0] == '"'  or item[0] == "'":  item = item[1:]
                if item[-1] == '"' or item[-1] == "'": item = item[:-1]
                items_.append(item)
            except:
                if item[0] == '"'  or item[0] == "'":  item = item[1:]
                if item[-1] == '"' or item[-1] == "'": item = item[:-1]
                items_.append('"'+item+'"')
        output_l.append(":".join(items_)) 
    print(", ".join(output_l))   
