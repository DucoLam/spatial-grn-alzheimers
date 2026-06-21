import fitz
d = fitz.open('/tudelft.net/staff-umbrella/ScReNI/bsc-screni/docs/paper.pdf')
out = []
for p in d:
    out.append(p.get_text())
open('/tudelft.net/staff-umbrella/ScReNI/dflam/_paper.txt','w').write('\n'.join(out))
print('FITZ pages', d.page_count)
