#=====================================file header=======================================#
# Author: Yifan Lin
# Contact: yil021@ucsd.edu
# Last update: 04/02/2021
# Purpose: to implement a user interactive time-series RNA-seq visualization toolkit with
#          Python Dash
#=======================================================================================#

import dash
import dash_bio as dashbio
import dash_html_components as html
import dash_core_components as dcc
import dash_table
import pandas as pd
import plotly.express as px
from dash.dependencies import Input, Output
import numpy
import base64
import datetime
import io

# define const: external stylesheet, groups, time points, colors, and string const
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
groups = ['all genes', 'differentially expressed genes']
tp = ['1hr', '2hr', '3hr', '4hr' , '6hr' , '9hr', '12hr', '16hr', '20hr', '24hr']
colors = {
    'background': '#ffffff',
    'text': '#111111'
}
clustergram_colors = ['#06081c', '#f49f76', '#e22e39', '#6a1e56', '#ffffff']
condition = 'OSvsPS'
p = '_Padjusted'
l = '_logFC'

# init variables
df = pd.read_csv('deseq2_all_genes.csv')
df_index = df.set_index('GeneName') # keep another indexed df for clustergram row names
columns = list(df_index.columns.values)
LFC = df_index[[i for i in columns if 'logFC' in i]]
df_pca = df[df['GeneName'].isin(['KLF4', 'KLF2', 'HR', 'NPR1', 'FGFR3', 'CXCR4'])] #p0.1; lfc1
tp_num_default = [1, 9, 24]

# func: DE_filter
# @param df: dataframe
# @param lfc: log2 FC cutoff
# @param padj: adjusted pval cutoff
# @param tpList: timepoint(s) the conditions applied to
def DE_filter(df, lfc, padj, tpList):
    for i in tp:
        tp_num = i[:-2]
        padjname = condition+tp_num+p
        lfcname = condition+tp_num+l
        df[str(tp_num)+'lfcpass'] = numpy.where(abs(df[lfcname]) > lfc, True, False)
        df[str(tp_num)+'padjpass'] = numpy.where(df[padjname] < padj, True, False)
    
    if tpList == ['all']:
        tpList = tp
    for i in tpList:
        df = df[df[i[:-2]+'lfcpass']==True]
        df = df[df[i[:-2]+'padjpass']==True]
    df = df.reset_index()
    del df['index']
    return df

# func: parse_contents
# @param contents: contents of upload (type and string)
# @param filename: name of file to parse
# @param date: last update date of file
def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
            df_index = df.set_index('GeneName') # keep another indexed df for clustergram row names
            columns = list(df_index.columns.values)
            LFC = df_index[[i for i in columns if 'logFC' in i]]
            df_pca = df[df['GeneName'].isin(['KLF4', 'KLF2', 'HR', 'NPR1', 'FGFR3', 'CXCR4'])] #p0.1; lfc1
            tp_num_default = [1, 9, 24]
            # partial display as there are too many rows
            df = df[0:25]
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
            
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

    return html.Div([
        html.I(filename),
        html.I(' preview'),
        html.Br(),
        html.I(datetime.datetime.fromtimestamp(date)),

        dash_table.DataTable(
            data=df.to_dict('records'),
            columns=[{'name': i, 'id': i} for i in df.columns]
        ),

        html.Hr(),  # horizontal line
    ])

# app major layout style and components
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.layout = html.Div(style={'backgroundColor': colors['background'], 'width': '70%','padding-left':'15%', 'padding-right':'15%'}, children=[ html.Br(),
    html.H1(
        children = 'Time-series RNA-seq Visualization Toolkit',
        style={
            'textAlign': 'center',
            'color': colors['text'],
            'fontSize': 45}),
    html.I('ver 1.0.0, 04/02/2021'), html.Br(), 
    html.I('ver 2.0.0, 04/04/2021'), html.Br(), html.Br(),
#========================pt0: upload a csv file to get started==========================#
    html.H1(
        children = '0 - Please upload your time-series RNA-seq dataset to get started',
        style={
            'textAlign': 'center',
            'color': colors['text'],
            'fontSize': 35}),
    dcc.Upload(
        id='upload-file',
        children=html.Div([html.A('Upload a File in .csv Format')]),
        filename='deseq2_all_genes.csv',
        style={
            'width': '100%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '10px'},
        # Allow multiple files to be uploaded
        multiple=True),
    html.Div(id='output-file-upload'),

#===========pt1: volcano plot to explore genes of interest based on custom LFC==========#
    html.H1(
        children = '1 - Volcanoplots of Transcriptional Levels at Inidividual Timepoint',
        style={
            'textAlign': 'center',
            'color': colors['text'],
            'fontSize': 35}),html.Br(),

    # dropdown to customize tp of volcano plots
   dcc.Dropdown(
        id='volcano-tp-input',
        options=[
            {'label': row, 'value': row} for row in tp],
        value='1hr'), html.Br(),
    
    html.I('Use the slider to customize Log2 Fold Change range: '),
    # slilder to customize LFC TODO: remove outlier and display
    dcc.RangeSlider(
        id='volcanoplot-input-lfc',
        min=-3,
        max=3,
        step=0.1,
        marks={
            i: {'label': str(i)} for i in range(-3, 3)
        },
        value=[-1.3, 1.3]), 
    html.I('Note: significant p-adj cutoff == 0.05 in this section. '), html.Hr(),

    # volcano plot div
    dcc.Graph(
        id='volcanoplot',
        figure=dashbio.VolcanoPlot(
            dataframe=df,
            snp="Entrez", # Entrez IDs
            gene="GeneName", # Gene Names
            p="OSvsPS6_Padjusted", # p-adj value column name
            genomewideline_value=0.05, # default p-adj cutoff
            genomewideline_color='#EF553B',
            effect_size="OSvsPS6_logFC", # effectSize column name (Log2FC)
            effect_size_line=[-1.3, 1.3], # default LFC cutoff
            effect_size_line_color='#AB63FA', 
            highlight_color='#34c6eb',
            col='#111111',
            title="Volcano plot of Human Umbilical Vein Cell RNAseq at 1hr",
            xlabel="log2 Fold Change",
            ylabel="-log10(p_adjusted)")), html.Hr(),html.Br(),

#=====================pt2: heatmap overview of transcriptional levels===================#
    html.H1(
        children = '2 - Transcriptional Levels Examination with Clustergarm',
        style={
            'textAlign': 'center',
            'color': colors['text'],
            'fontSize': 35}),
    html.I(
        children = 'Guide:',
        style={'fontSize': 20}), html.Br(),
    # user guide
    html.I('1) Choose a subset of genes you wish to see with the dropdown'), html.Br(),
    html.I('2) Type your log2 fold change and p-adjusted value cutoffs to customize your search, and click apply'), html.Br(),
    html.I('3) Specify which timepoint(s) of transcriptional data you wish to examine with checkboxes; this section shows genes applied to above cutoffs in all of your chosen timepoints. '), html.Br(), html.Br(),
    html.I("Please choose a subset of genes of which mRNA levels you would like to see. "), html.Br(),
    html.I("Note: due to limited memory, max of 1000 genes are to display in this page, and please wait for loading patiently. "), html.Br(),
    
    # clustergram dropdown options
    dcc.Dropdown(
        id='clustergram-input',
        options=[
            {'label': group, 'value': group} for group in groups],
        value=groups[0], 
        style={
            "width":"100%", 
            "align-items": "center", 
            "justify-content": "center"}), html.Br(), html.Br(),
    html.I("The following cutoffs only apply when you decide to see differentially expressed genes. "), html.Br(),

    # cutoffs input and submit
    html.I("Please input your log2 fold change cutoff (step = 0.1): "),
    dcc.Input(
        id = 'lfcinput', 
        type = 'number', 
        value = 1.3,
        step = 0.1),
    html.Button(
        'apply',
        id = 'applyLFC',
        n_clicks = 0
        ), html.Br(),
    html.I("Please input your p-adjusted value cutoff (step = 0.01): "),
    dcc.Input(
        id = 'padjinput', 
        type = 'number', 
        value = 0.05,
        step = 0.01),
    html.Button(
        'apply',
        id = 'applyPadj',
        n_clicks = 0
        ), html.Br(),

    # timepoint checkbox
    html.I("Above cutoffs apply to timepoint(s): "),
    html.Div(
        dcc.Checklist(
            id='clustergram-tp-input',
            options=[
                {'label': '1hr', 'value': '1hr'},
                {'label': '2hr', 'value': '2hr'},
                {'label': '3hr', 'value': '3hr'},
                {'label': '4hr', 'value': '4hr'},
                {'label': '6hr', 'value': '6hr'},
                {'label': '9hr', 'value': '9hr'},
                {'label': '12hr', 'value': '12hr'},
                {'label': '16hr', 'value': '16hr'},
                {'label': '20hr', 'value': '20hr'},
                {'label': '24hr', 'value': '24hr'},
                {'label': 'all', 'value': 'all'}],         
            value=[],
            labelStyle={'display': 'inline-block'}),
        style={'display': 'inline-block'}), 
    # clustergram div
    dcc.Graph(
        id='my-clustergram',
        figure=dashbio.Clustergram(
            row_labels=list(df['GeneName'][0:1000]),
            column_labels=[i for i in columns if 'logFC' in i],
            data=LFC.loc[list(df['GeneName'][0:1000])].values,
            color_threshold={
                'row': -10,
                'col': 10},
            height=1200,
            width=900,
            color_map= [[0.0, '#06081c'],[0.25, '#f49f76'],[0.5, '#e22e39'],[0.75, '#6a1e56'],[1, '#ffffff']])), 
    html.Hr(), html.Br(),

#=====================pt3: PCA analysis based on top ranked DE genes====================#
    html.H1(
        children = '3 - Customized 3D Principle Component Analysis',
        style={
            'textAlign': 'center',
            'color': colors['text'],
            'fontSize': 35}), html.Br(),

    # customize 3 timepoints in pca plots
    html.I("Please choose 3 timepoints as dimensions of the analysis: "),
    dcc.Checklist(
        id='pca-tp-input',
        options=[
            {'label': '1hr', 'value': '1hr'},
            {'label': '2hr', 'value': '2hr'},
            {'label': '3hr', 'value': '3hr'},
            {'label': '4hr', 'value': '4hr'},
            {'label': '6hr', 'value': '6hr'},
            {'label': '9hr', 'value': '9hr'},
            {'label': '12hr', 'value': '12hr'},
            {'label': '16hr', 'value': '16hr'},
            {'label': '20hr', 'value': '20hr'},
            {'label': '24hr', 'value': '24hr'}
        ],
        value=['1hr', '9hr', '24hr']), html.Br(),
    
    # customize genes displayed in pca plot
    html.I("Please choose gene(s) of your interests to visualize in PCA plot: "),
    dcc.Dropdown(
            id='pca-gene-input',
            options=[
                {'label': i, 'value': i} for i in list(df_index.index)],
            value=['KLF4', 'KLF2', 'HR', 'NPR1', 'FGFR3', 'CXCR4'],
            multi=True), 
    # 3D pca div
    dcc.Graph(
        id='pca-plot',
        figure= px.scatter_3d(
            df_pca, 
            x='OSvsPS'+str(tp_num_default[0])+'_logFC', 
            y='OSvsPS'+str(tp_num_default[1])+'_logFC', 
            z='OSvsPS'+str(tp_num_default[2])+'_logFC',
            labels={
                'OSvsPS'+str(tp_num_default[0])+'_logFC':'1hr log2FC', 
                'OSvsPS'+str(tp_num_default[1])+'_logFC':'9hr log2FC', 
                'OSvsPS'+str(tp_num_default[2])+'_logFC':'24hr log2FC'},
            color='GeneName',
            width=800,
            height=600)), html.Hr()
])

# file upload callback and update function
# @param contents: upload info including type and string
# @param filenames: filename
# @param lastMods: lastly modified time
@app.callback(
    dash.dependencies.Output('output-file-upload', 'children'),
    dash.dependencies.Input('upload-file', 'contents'),
    dash.dependencies.State('upload-file', 'filename'),
    dash.dependencies.State('upload-file', 'last_modified'))
def upload_file(contents, filenames, lastMods):
    if contents is not None:
        children = [
            parse_contents(c, n, d) for c, n, d in
            zip(contents, filenames, lastMods)]
        return children

# volcano plot callback and update function
# @param: 'volcano-tp-input' tp
# @param: 'volcanoplot-input-lfc' lfc
# @return: 'volcanoplot' 
@app.callback(
    dash.dependencies.Output('volcanoplot', 'figure'),
    [dash.dependencies.Input('volcano-tp-input', 'value'),
    dash.dependencies.Input('volcanoplot-input-lfc', 'value')])
def update_volcanoplot(tp, lfc):
    if tp == None:
        tp = '1hr'
    padjCol = 'OSvsPS'+tp[0]+'_Padjusted'
    lfcCol = 'OSvsPS'+tp[0]+'_logFC'
    tp_num = tp[:-2]
        
    # print((lfcCol))
    return dashbio.VolcanoPlot(
        dataframe=df,
        #dataframe=df[(df[tp_num+"pass"] == True)],
        snp="Entrez", # Entrez IDs
        gene="GeneName", # Gene Names
        p=padjCol, # p-adj value column name
        genomewideline_value=1.3, # default p-adj cutoff
        genomewideline_color='#EF553B',
        effect_size=lfcCol, # effectSize column name (Log2FC)
        effect_size_line=lfc, # updated effect size
        effect_size_line_color='#AB63FA', 
        highlight_color='#34c6eb',
        col='#111111',
        title="Volcano plot of Human Umbilical Vein Cell RNAseq at "+tp[0]+"hr",
        xlabel="log2 Fold Change",
        ylabel="-log10(p_adjusted)"
    )

# clustergram callback and update function
# @param: 'clustergram-input' group
# @param: 'applyLFC' lfcClicks
# @param: 'applyPadj' padjClicks
# @param: 'clustergram-tp-input' tp
# @param: 'lfcinput' lfc
# @param: 'padjinput' padj
# @return: 'my-clustergram'
@app.callback(
    dash.dependencies.Output('my-clustergram', 'figure'),
    [dash.dependencies.Input('clustergram-input', 'value')],
    [dash.dependencies.Input('applyLFC', 'n_clicks')],
    [dash.dependencies.Input('applyPadj', 'n_clicks')],
    [dash.dependencies.Input('clustergram-tp-input', 'value')],
    [dash.dependencies.State('lfcinput', 'value')],
    [dash.dependencies.State('padjinput', 'value')])
def update_clustergram(group, lfcClicks, padjClicks, tp, lfc, padj):    
    displaydf = list()
    if group == groups[0]:
        displaydf = list(df['GeneName'][0:1000])
    
    elif group == groups[1] and len(tp) == 0 :
        return {
        "layout": {
            "xaxis": {
                "visible": False},
            "yaxis": {
                "visible": False},
            "annotations": [{
                    "text": "Please select a subgroup or at least a timepoint to visualize data",
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {
                        "size": 15}}]}}

    
    elif group == groups[1] and len(tp) != 0 :
        newdf = DE_filter(df, lfc, padj, tp)
        displaydf = list(newdf['GeneName'])

    elif 'all' in tp:
        tp_temp = ['1hr', '2hr', '3hr', '4hr' , '6hr' , '9hr', '12hr', '16hr', '20hr', '24hr']
        newdf = DE_filter(df, lfc, padj, tp_temp)
        displaydf = list(newdf['GeneName'])

    fig = dashbio.Clustergram(
        row_labels=displaydf,
        column_labels=[i for i in columns if l in i],
        data=LFC.loc[displaydf].values,
        color_threshold={
            'row': -10,
            'col': 10
        },
        height=1200,
        width=900,
        color_map= [
        [0.0, clustergram_colors[0]],
        [0.25, clustergram_colors[1]],
        [0.5, clustergram_colors[2]],
        [0.75, clustergram_colors[3]],
        [1, clustergram_colors[4]]
        ])
    return fig

# pca callback and update function
# @param: 'pca-tp-input' tp
# @param: 'pca-gene-input' genes
# @return: 'pca-plot'
@app.callback(
    dash.dependencies.Output('pca-plot', 'figure'),
    [dash.dependencies.Input('pca-tp-input', 'value')],
    [dash.dependencies.Input('pca-gene-input', 'value')])
def update_pca(tp, genes):
    if len(tp) != 3: # return an empty graph
        return {
        "layout": {
            "xaxis": {
                "visible": False},
            "yaxis": {
                "visible": False},
            "annotations": [{
                    "text": "Please select exactly 3 timepoints",
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {
                        "size": 15}}]}}
    else:
        tp_num=[]
        for i in tp:
            tp_num.append(int(i[:-2]))

    df_pca = df[df['GeneName'].isin(genes)]
    figure= px.scatter_3d(
        df_pca, 
        x=condition+str(tp_num[0])+l, 
        y=condition+str(tp_num[1])+l, 
        z=condition+str(tp_num[2])+l,
        labels={
            condition+str(tp_num[0])+l:str(tp_num[0])+'hr log2FC', 
            condition+str(tp_num[1])+l:str(tp_num[1])+'hr log2FC', 
            condition+str(tp_num[2])+l:str(tp_num[2])+'hr log2FC'},
        color='GeneName',
        width=800,
        height=600)
    return figure

if __name__ == '__main__':
    app.run_server(debug=True)