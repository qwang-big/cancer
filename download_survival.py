import requests, urllib, json, sys
d={'methodoption': 'os', 
	'dataset': 'LUAD', 
	'signature': 'kras', 
	'highcol': '%23ff0000', 
	'lowcol': '%230000ff', 
	'groupcutoff1': '50', 
	'groupcutoff2': '50', 
	'axisunit': 'month', 
	'ifhr': 'hr', 
	'ifconf': 'conf', 
	'signature_norm': '', 
	'is_sub': 'false', 
	'subtype': ''}
s=sys.argv[1]
d['signature']=s
r=requests.post("http://gepia2.cancer-pku.cn/assets/PHP4/survival_zf.php", d)
urllib.request.urlretrieve('http://gepia2.cancer-pku.cn/tmp/'+json.loads(r.text[2:])['outdir'],s+'.pdf')
