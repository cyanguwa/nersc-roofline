import pandas as pd
file='my-vtune.knl.summary'
df=pd.read_csv(file)

total_read=df.filter(regex='UNC_M_CAS_COUNT.RD').sum(axis=1)*64
total_write=df.filter(regex='UNC_M_CAS_COUNT.WR').sum(axis=1)*64
print('--->DDR Report')
print('--->Total Bytes read = '+str(total_read[0]))
print('--->Total Bytes written = '+str(total_write[0]))
print('--->Total Bytes = '+str(total_read[0] + total_write[0] ))

total_read=df.filter(regex='UNC_E_RPQ_INSERTS').sum(axis=1)*64
total_write=df.filter(regex='UNC_E_WPQ_INSERTS').sum(axis=1)*64
print('--->MCDRAM Report')
print('--->Total Bytes read = '+str(total_read[0]))
print('--->Total Bytes written = '+str(total_write[0]))
print('--->Total Bytes = '+str(total_read[0] + total_write[0]))


