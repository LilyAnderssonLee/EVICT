#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import base64
from io import BytesIO
import datetime
import glob
import os
import sys
import argparse
from jinja2 import Environment, BaseLoader

def main():
    parser = argparse.ArgumentParser(
        description='Generate Enterovirus typing reports from BLAST results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
        Examples:
        python EV_typing.py --ticket_nr 1003460 --blast-file /path/to/sample.blast
        """
    )
    parser.add_argument ('--ticket-nr', type=int, required=True, help='Ticket number')
    parser.add_argument('--blast-file', type=str, required=True, help='Path to the BLAST result file (.blast)')
    parser.add_argument('--output-dir', type=str, default='results', help='Base output directory')
    parser.add_argument('--dpi', type=int, default=300, help='DPI for output plots')

    args = parser.parse_args()
    ticket_nr, blast_file, output_base, plot_dpi = args.ticket_nr, args.blast_file, args.output_dir, args.dpi

    if not os.path.exists(blast_file):
        print(f"❌ Error: BLAST file not found: {blast_file}")
        sys.exit(1)

    # Extract file name (without extension) to use as sequence name
    file_name = os.path.basename(blast_file)
    seq_name = file_name.replace('.blast', '')

    # Ensure output dir exists
    os.makedirs(output_base, exist_ok=True)

    # Initialize Jinja2 environment
    jinja_env = Environment(loader=BaseLoader())

    # Unified Jinja2 template for both regular and error reports
    UNIFIED_TEMPLATE = """<!DOCTYPE html>
<html lang="sv">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Enterovirus Rapport - {{ seq_name }}</title>
    <style>
        :root {
            --primary-color: #2c3e50;
            --secondary-color: #3498db;
            --accent-color: #e74c3c;
            --success-color: #27ae60;
            --warning-color: #f39c12;
            --light-gray: #ecf0f1;
            --dark-gray: #7f8c8d;
            --text-color: #2c3e50;
            --border-color: #bdc3c7;
            --shadow: 0 2px 10px rgba(0,0,0,0.1);
        }

        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: var(--text-color);
            background-color: #f8f9fa;
            margin: 0;
            padding: 20px;
        }

        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 12px;
            box-shadow: var(--shadow);
            overflow: hidden;
        }

        .header {
            background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
            color: white;
            padding: 2rem;
            text-align: center;
            position: relative;
        }

        .header::before {
            content: '';
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background: url('data:image/svg+xml,<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 20"><defs><pattern id="grid" width="10" height="10" patternUnits="userSpaceOnUse"><path d="M 10 0 L 0 0 0 10" fill="none" stroke="white" stroke-width="0.5" opacity="0.1"/></pattern></defs><rect width="100" height="20" fill="url(%23grid)"/></svg>');
        }

        .header-content {
            position: relative;
            z-index: 1;
        }

        .header h1 {
            font-size: 2em;
            font-weight: 300;
            margin-bottom: 0.5rem;
            text-shadow: 0 2px 4px rgba(0,0,0,0.3);
        }

        .header p {
            font-size: 1.5em;
            font-weight: 400;
            margin: 0;
        }

        .metadata {
            background: var(--light-gray);
            padding: 1.5rem 2rem;
            border-bottom: 1px solid var(--border-color);
        }

        .metadata-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1rem;
        }

        .metadata-item {
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }

        .metadata-label {
            font-weight: 600;
            color: var(--dark-gray);
            min-width: 80px;
        }

        .metadata-value {
            font-weight: 500;
            color: var(--primary-color);
        }

        .content {
            padding: 2rem;
        }

        .description {
            background: #f8f9fb;
            border-left: 4px solid var(--secondary-color);
            padding: 1.5rem;
            margin-bottom: 2rem;
            border-radius: 0 8px 8px 0;
        }

        .warning {
            background: linear-gradient(135deg, #fff3cd, #ffeaa7);
            border: 1px solid var(--warning-color);
            border-radius: 8px;
            padding: 1rem;
            margin: 1.5rem 0;
            display: flex;
            align-items: center;
            gap: 1rem;
        }

        .warning-icon {
            font-size: 1.5em;
            color: var(--warning-color);
        }

        .warning-text {
            font-weight: 500;
            color: #856404;
        }

        .figure-container {
            background: white;
            border: 1px solid var(--border-color);
            border-radius: 12px;
            padding: 2rem;
            margin: 2rem 0;
            box-shadow: 0 4px 6px rgba(0,0,0,0.05);
        }

        .figure-container img {
            width: 100%;
            height: auto;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }

        .figure-legend {
            background: var(--light-gray);
            border-radius: 8px;
            padding: 1.5rem;
            margin-top: 1rem;
        }

        .figure-legend h3 {
            color: var(--primary-color);
            margin-bottom: 1rem;
            font-size: 1.1em;
            font-weight: 600;
        }

        .figure-legend ul {
            list-style: none;
            padding: 0;
        }

        .figure-legend li {
            padding: 0.5rem 0;
            border-bottom: 1px solid #ddd;
            position: relative;
            padding-left: 1.5rem;
        }

        .figure-legend li:last-child {
            border-bottom: none;
        }

        .figure-legend li::before {
            content: '📊';
            position: absolute;
            left: 0;
            top: 0.5rem;
        }

        .links-section {
            background: var(--light-gray);
            border-radius: 8px;
            padding: 1rem;
            margin: 1.5rem 0;
        }

        .links-section h3 {
            color: var(--primary-color);
            margin-bottom: 1rem;
            font-size: 1.1em;
            text-align: center;
        }

        .links-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 0.8rem;
        }

        .link-card {
            background: white;
            border-radius: 6px;
            padding: 0.8rem;
            text-decoration: none;
            color: var(--text-color);
            transition: all 0.3s ease;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            border-left: 3px solid var(--secondary-color);
        }

        .link-card:hover {
            transform: translateY(-1px);
            box-shadow: 0 2px 6px rgba(0,0,0,0.15);
            text-decoration: none;
            color: var(--secondary-color);
        }

        .link-card h4 {
            font-size: 0.9em;
            font-weight: 600;
            margin-bottom: 0.3rem;
        }

        .link-card p {
            font-size: 0.85em;
            color: var(--dark-gray);
            margin: 0;
        }

        .external-tool {
            background: var(--light-gray);
            color: var(--text-color);
            padding: 1rem;
            border-radius: 8px;
            margin: 1.5rem 0;
            text-align: center;
            border: 1px solid var(--border-color);
        }

        .external-tool a {
            color: var(--secondary-color);
            text-decoration: none;
            font-weight: 500;
            font-size: 1em;
        }

        .external-tool a:hover {
            text-decoration: underline;
        }

        .contigs-section {
            margin-top: 2rem;
        }

        .contigs-section h3 {
            color: var(--primary-color);
            margin-bottom: 1rem;
            font-size: 1.3em;
            padding-bottom: 0.5rem;
            border-bottom: 2px solid var(--secondary-color);
        }

        .contigs-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            background: var(--light-gray);
            padding: 1rem 1.5rem;
            border-radius: 8px 8px 0 0;
            border-bottom: 2px solid var(--border-color);
        }

        .contigs-summary {
            background: #2c3e50;
            color: #ecf0f1;
            padding: 1.5rem;
            border-radius: 0;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
            line-height: 1.4;
            overflow-x: auto;
            box-shadow: inset 0 2px 4px rgba(0,0,0,0.2);
        }

        .contigs-full {
            background: #2c3e50;
            color: #ecf0f1;
            padding: 1.5rem;
            border-radius: 0 0 8px 8px;
            font-family: 'Courier New', monospace;
            font-size: 0.9em;
            line-height: 1.4;
            overflow-x: auto;
            box-shadow: inset 0 2px 4px rgba(0,0,0,0.2);
            display: none;
        }

        .contigs-container {
            border: 1px solid var(--border-color);
            border-radius: 8px;
            overflow: hidden;
        }

        .toggle-btn {
            background: var(--secondary-color);
            color: white;
            border: none;
            padding: 0.5rem 1rem;
            border-radius: 6px;
            cursor: pointer;
            font-size: 0.9em;
            font-weight: 500;
            transition: all 0.3s ease;
            display: flex;
            align-items: center;
            gap: 0.5rem;
        }

        .toggle-btn:hover {
            background: #2980b9;
            transform: translateY(-1px);
        }

        .toggle-btn:active {
            transform: translateY(0);
        }

        .toggle-icon {
            transition: transform 0.3s ease;
        }

        .toggle-btn.expanded .toggle-icon {
            transform: rotate(180deg);
        }

        .footer {
            background: var(--primary-color);
            color: white;
            text-align: center;
            padding: 1rem;
            font-size: 0.9em;
        }

        @media (max-width: 768px) {
            body {
                padding: 10px;
            }

            .header h1 {
                font-size: 1.8em;
            }

            .content {
                padding: 1rem;
            }

            .metadata-grid {
                grid-template-columns: 1fr;
            }

            .links-grid {
                grid-template-columns: 1fr;
            }

            .contigs-header {
                flex-direction: column;
                gap: 1rem;
                align-items: stretch;
            }

            .toggle-btn {
                justify-content: center;
            }
        }

        @media print {
            body {
                background: white;
                padding: 0;
            }

            .container {
                box-shadow: none;
                border-radius: 0;
            }

            .link-card:hover {
                transform: none;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }

            .toggle-btn {
                display: none;
            }

            .contigs-full {
                display: block !important;
            }
        }
    </style>
    <script>
        function toggleContigs() {
            const summaryDiv = document.getElementById('contigs-summary');
            const fullDiv = document.getElementById('contigs-full');
            const toggleBtn = document.getElementById('toggle-btn');
            const toggleIcon = document.querySelector('.toggle-icon');

            if (fullDiv.style.display === 'none' || fullDiv.style.display === '') {
                summaryDiv.style.display = 'none';
                fullDiv.style.display = 'block';
                toggleBtn.innerHTML = '<span class="toggle-icon">▲</span> Dölj innehåll';
                toggleBtn.classList.add('expanded');
            } else {
                summaryDiv.style.display = 'block';
                fullDiv.style.display = 'none';
                toggleBtn.innerHTML = '<span class="toggle-icon">▼</span> Visa innehåll';
                toggleBtn.classList.remove('expanded');
            }
        }
    </script>
</head>
<body>
    <div class="container">
        <header class="header">
            <div class="header-content">
                <h1>Enterovirus Genotypning</h1>
                <p>Prov: {{ seq_name }}</p>
            </div>
        </header>

        <div class="metadata">
            <div class="metadata-grid">
                <div class="metadata-item">
                    <span class="metadata-label">Ticket:</span>
                    <span class="metadata-value">{{ ticket_nr }}</span>
                </div>
                <div class="metadata-item">
                    <span class="metadata-label">Utförd:</span>
                    <span class="metadata-value">{{ time_stamp }}</span>
                </div>
            </div>
        </div>

        <div class="content">
            <div class="description">
                {% if is_error_report %}
                <p><strong>Sammanställning av resultat</strong> för alla contigs (ofiltrerade data).</p>
                {% else %}
                <p><strong>Sammanställning av resultat</strong> för contigs med täckning minst 50x och längd minst 200bp.</p>
                {% endif %}
            </div>

            {% if is_error_report %}
            <div class="warning">
                <span class="warning-icon">⚠️</span>
                <span class="warning-text">OBS! Inga contigs uppfyllde kriterierna för standardanalys (minst 200bp längd och 50x täckning). Visar resultat för alla tillgängliga contigs.</span>
            </div>
            {% endif %}

            {% if varningstext and not is_error_report %}
            <div class="warning">
                <span class="warning-icon">⚠️</span>
                <span class="warning-text">{{ varningstext | safe }}</span>
            </div>
            {% endif %}

            <div class="figure-container">
                <img src="{{ img_data_uri }}" alt="BLAST analys resultat">

                <div class="figure-legend">
                    <h3>Förklaring av figurer</h3>
                    <ul>
                        <li><strong>Figur 1 och 2:</strong> Cirklar anger identitet per matchad sekvens, romber anger medelvärde per contig.</li>
                        <li><strong>Identitet:</strong> Hur många nukleotider i alignment är identiska.</li>
                        <li><strong>Bit score:</strong> Normaliserad score för att statistiskt redovisa hur bra alignment är.</li>
                    </ul>
                </div>
            </div>

            <div class="links-section">
                <h3>Referensinformation - Enterovirus genotyper per species</h3>
                <div class="links-grid">
                    <a href="https://www.picornaviridae.com/ensavirinae/enterovirus/ev-a/ev-a.htm" class="link-card">
                        <h4>Enterovirus Species A</h4>
                        <p>Komplett lista över genotyper</p>
                    </a>
                    <a href="https://www.picornaviridae.com/ensavirinae/enterovirus/ev-b/ev-b.htm" class="link-card">
                        <h4>Enterovirus Species B</h4>
                        <p>Komplett lista över genotyper</p>
                    </a>
                    <a href="https://www.picornaviridae.com/ensavirinae/enterovirus/ev-c/ev-c.htm" class="link-card">
                        <h4>Enterovirus Species C</h4>
                        <p>Komplett lista över genotyper</p>
                    </a>
                    <a href="https://www.picornaviridae.com/ensavirinae/enterovirus/ev-d/ev-d.htm" class="link-card">
                        <h4>Enterovirus Species D</h4>
                        <p>Komplett lista över genotyper</p>
                    </a>
                </div>
            </div>

            <div class="external-tool">
                <p>För ytterligare analys kan contigerna analyseras med:</p>
                <a href="https://mpf.rivm.nl/mpf/typingtool/enterovirus/">RIVM Enterovirus Genotyping Tool</a>
            </div>

            <div class="contigs-section">
                <h3>Contiger</h3>
                <div class="contigs-container">
                    <div class="contigs-header">
                        <div>
                            <strong>🧬 Sekvensdata:</strong>
                            <span style="color: var(--dark-gray); margin-left: 1rem;">{{ contig_count }} contig(er) funna</span>
                        </div>
                        <button class="toggle-btn" id="toggle-btn" onclick="toggleContigs()">
                            <span class="toggle-icon">▼</span> Visa innehåll
                        </button>
                    </div>
                    <div class="contigs-summary" id="contigs-summary">{{ contigs_summary | safe }}</div>
                    <div class="contigs-full" id="contigs-full">{{ contigs_content | safe }}</div>
                </div>
            </div>
        </div>

        <footer class="footer">
            <p>Genererad av Enterovirus Genotypning System: Klinisk Mikrobiologi, Karolinska Universitetssjukhuset</p>
        </footer>
    </div>
</body>
</html>"""

    def extract_contig_headers(fasta_content):
        """Extract only the header lines from FASTA content"""
        headers = []
        for line in fasta_content.split('\n'):
            if line.startswith('>'):
                headers.append(line)
        return headers

    def generate_normal_report_data(file_path):
        """Generate data for normal reports (filtered)"""
        file_name = os.path.basename(file_path)
        seq_name = file_name.replace('.blast', '')

        # Read file
        df = pd.read_csv(file_path,
                         header=0)

        # Split qseqid to obtain contig name and coverage
        df[['contig','temp1']] = df['qseqid'].str.split('_length_',expand=True)
        df[['length','coverage']] = df['temp1'].str.split('_cov_',expand=True)

        # Remove temporary and redundant columns
        df = df.drop(['qseqid', 'temp1', 'length'], axis=1)

        #convert coverage to numeric
        df["coverage"] = pd.to_numeric(df["coverage"])

        # Remove contigs with less than 200bp length and 50x coverage
        df = df.loc[df['qlen'] > 200]
        df = df.loc[df['coverage'] > 50]

        # CHECK: If no contigs remain after filtering, raise an exception
        if df.empty or df['contig'].nunique() == 0:
            raise ValueError(f"No contigs meet filtering criteria (≥200bp length, ≥50x coverage) for {file_name}")

        # Check if filtered FASTA file exists
        seq_file = f"{output_base}/{ticket_nr}/ev_contig/{file_name.replace('.blast', '_200bp_minCov50.fasta')}"
        if not os.path.exists(seq_file):
            raise FileNotFoundError(f"Filtered FASTA file not found: {seq_file}")

        # Count number of contigs
        n_contigs = df['contig'].nunique()

        # Coverage of contigs
        u_contigs = df[['contig','coverage','scomname']].drop_duplicates()

        # Count number of genotypes
        n_scomname = df['scomname'].nunique()

        # Length of contigs
        l_contigs = df[['contig','qlen','scomname']].drop_duplicates()

        # Additional safety checks
        if n_contigs == 0 or n_scomname == 0:
            raise ValueError(f"Invalid contig or genotype count for {file_name}")

        pident_means = df.groupby('scomname')['pident'].mean().sort_values(ascending=True)
        bitscore_means = df.groupby('scomname')['bitscore'].mean().sort_values(ascending=True)
        coverage_means = df.groupby('contig')['coverage'].mean().sort_values(ascending=True)
        length_means = df.groupby('contig')['qlen'].mean().sort_values(ascending=True)

        pident_order = pident_means.index.tolist()
        bitscore_order = bitscore_means.index.tolist()
        coverage_order = coverage_means.index.tolist()
        length_order = length_means.index.tolist()

        # Set plot parameters
        if n_scomname+n_contigs == 2:
            fig_height = 3
        elif n_scomname+n_contigs > 8:
            fig_height = (n_scomname+n_contigs)*0.8
        else:
            fig_height = n_scomname+n_contigs

        if n_contigs == 1:
            h_ratio = 1
        else:
            h_ratio = n_contigs * 0.7

        if n_contigs > 1:
            legend_status = True
            fig_width = 14
        else:
            legend_status = False
            fig_width = 10

        fig = plt.figure(figsize=(fig_width, fig_height))
        gs = fig.add_gridspec(2, 2, height_ratios=[n_scomname, h_ratio])

        # Plot 1: identity score
        ax = fig.add_subplot(gs[0 , 0])
        ax.set_title("Figur 1: Identitet per genotyp")

        # Define contents of plot
        x_cat = "pident"
        y_cat = "scomname"
        hue_cat = "contig"

        sns.pointplot(
            data=df,
            x=x_cat,
            y=y_cat,
            hue=hue_cat,
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", zorder=10, order = pident_order, legend =False
        )

        sns.stripplot(
            data=df,
            x=x_cat,
            y=y_cat,
            hue=hue_cat,
            dodge=True, alpha=.4, legend=False, jitter = 0.3, order = pident_order
        )

        if n_scomname > 1:
            for i in range(n_scomname):
                if i % 2 == 1:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)

        ax.set_ylim(-0.5, n_scomname-0.3)
        ax.grid(True, axis="x", linestyle="--")
        ax.set(xlabel='BLAST identitet (%)', ylabel='Genotyp')

        # Plot 2: bit score
        ax = fig.add_subplot(gs[0 , 1])
        ax.set_title("Figur 2: Bit score per genotyp")

        x_cat = "bitscore"
        y_cat = "scomname"
        hue_cat = "contig"

        sns.pointplot(
            data=df,
            x=x_cat,
            y=y_cat,
            hue=hue_cat,
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", zorder=10, order = bitscore_order, legend = legend_status
        )

        sns.stripplot(
            data=df,
            x=x_cat,
            y=y_cat,
            hue=hue_cat,
            dodge=True, alpha=.4, legend=False, jitter = 0.3, order = bitscore_order
        )
        if n_scomname > 1:
            for i in range(n_scomname):
                if i % 2 == 1:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)

        ax.set_ylim(-0.5, n_scomname-0.3)
        ax.grid(True, axis="x", linestyle="--")

        if legend_status == True:
            ax.legend().set_title("Contig")
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        ax.set(xlabel='BLAST bit score', ylabel=' ')

        # Plot 3: coverage
        ax = fig.add_subplot(gs[1, 0])
        ax.set_title("Figur 3: Täckning per contig")

        x_cat = "coverage"
        y_cat = "contig"
        hue_cat = y_cat

        sns.pointplot(
            data=u_contigs,
            x=x_cat,
            y=y_cat,
            hue=hue_cat,
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", order = coverage_order
        )

        if n_contigs > 1:
            for i in range(n_contigs):
                if i % 2 == 0:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)

        ax.set_ylim(-0.5, n_contigs-0.3)
        ax.grid(True, axis="x", linestyle="--")
        ax.set(xlabel='Täckning (x)', ylabel='Contig')
        ax.set_xscale('log')
        ax.set_xlim(left=50, right = 10*df['coverage'].max() )

        # Plot 4 length
        ax = fig.add_subplot(gs[1, 1])
        ax.set_title("Figur 4: Längd per contig")

        x_cat = "qlen"
        y_cat = "contig"
        hue_cat = y_cat

        sns.pointplot(
            data=l_contigs,
            x=x_cat,
            y=y_cat,
            hue=hue_cat,
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", order = length_order, legend = False
        )

        if n_contigs > 1:
            for i in range(n_contigs):
                if i % 2 == 0:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)

        ax.set_ylim(-0.5, n_contigs-0.3)
        ax.grid(True, axis="x", linestyle="--")
        ax.set(xlabel='Längd (bp)', ylabel=' ')
        ax.set_xlim(left=200, right = 1.1*df['qlen'].max())

        fig.subplots_adjust(hspace = 20, wspace=20)
        fig.tight_layout()

        # Convert plot to base64 with high resolution
        img_buffer = BytesIO()
        fig.savefig(img_buffer, format='png', dpi=plot_dpi, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        img_buffer.seek(0)
        img_base64 = base64.b64encode(img_buffer.read()).decode('utf-8')
        img_data_uri = f"data:image/png;base64,{img_base64}"

        # Read fasta-file (filtered)
        seq_file = f"{output_base}/{ticket_nr}/ev_contig/{file_name.replace('.blast', '_200bp_minCov50.fasta')}"
        with open(seq_file, 'r') as contigs_file:
            contigs_content = contigs_file.read()

        # Extract contig headers for summary
        contig_headers = extract_contig_headers(contigs_content)
        contig_count = len(contig_headers)

        # Create summary (just headers)
        contigs_summary = '<br>'.join(contig_headers)

        # Convert full content newlines to HTML breaks for display
        contigs_content_html = contigs_content.replace('\n', '<br>')

        # Warning for low identity
        varningstext = ''
        if df['pident'].max() < 99:
            varningstext = "&#9888; OBS! Högsta identitet har ett lågt värde (identitet < 99%) - fördjupad utredning nödvändig."

        return img_data_uri, contigs_content_html, contigs_summary, contig_count, varningstext

    def generate_error_report_data(file_path):
        """Generate data for error reports (unfiltered)"""
        file_name = os.path.basename(file_path)
        seq_name = file_name.replace('.blast', '')

        # Initialize default values in case of failure
        img_data_uri = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=="
        contigs_content = ""
        contigs_summary = ""
        contig_count = 0

        try:
            # Read BLAST file
            df = pd.read_table(file_path,
                               sep=",",
                              header=0)

            # Split qseqid to obtain contig name and coverage
            df[['contig','temp1']] = df['qseqid'].str.split('_length_',expand=True)
            df[['length','coverage']] = df['temp1'].str.split('_cov_',expand=True)

            # Remove temporary and redundant columns
            df = df.drop(['qseqid', 'temp1', 'length'], axis=1)

            # Convert coverage to numeric
            df["coverage"] = pd.to_numeric(df["coverage"])

            # Check if dataframe is empty or has no valid data
            if df.empty or df['contig'].nunique() == 0:
                print(f"Warning: No valid data found in BLAST file {file_name}")
                return img_data_uri, contigs_content, contigs_summary, contig_count

            # DON'T FILTER - use all contigs for error reports
            # Count number of contigs
            n_contigs = df['contig'].nunique()

            # Coverage of contigs
            u_contigs = df[['contig','coverage','scomname']].drop_duplicates()

            # Count number of genotypes
            n_scomname = df['scomname'].nunique()

            # Length of contigs
            l_contigs = df[['contig','qlen','scomname']].drop_duplicates()

            # Additional safety checks
            if n_contigs == 0 or n_scomname == 0:
                print(f"Warning: Invalid contig or genotype count for {file_name}")
                return img_data_uri, contigs_content, contigs_summary, contig_count

            # Generate the same plots as in the main function
            pident_means = df.groupby('scomname')['pident'].mean().sort_values(ascending=True)
            bitscore_means = df.groupby('scomname')['bitscore'].mean().sort_values(ascending=True)
            coverage_means = df.groupby('contig')['coverage'].mean().sort_values(ascending=True)
            length_means = df.groupby('contig')['qlen'].mean().sort_values(ascending=True)

            pident_order = pident_means.index.tolist()
            bitscore_order = bitscore_means.index.tolist()
            coverage_order = coverage_means.index.tolist()
            length_order = length_means.index.tolist()

            # Set plot parameters (same logic as main function)
            if n_scomname+n_contigs == 2:
                fig_height = 3
            elif n_scomname+n_contigs > 8:
                fig_height = (n_scomname+n_contigs)*0.8
            else:
                fig_height = n_scomname+n_contigs

            if n_contigs == 1:
                h_ratio = 1
            else:
                h_ratio = n_contigs * 0.7

            if n_contigs > 1:
                legend_status = True
                fig_width = 14
            else:
                legend_status = False
                fig_width = 10

            fig = plt.figure(figsize=(fig_width, fig_height))
            gs = fig.add_gridspec(2, 2, height_ratios=[n_scomname, h_ratio])

            # Plot 1: identity score
            ax = fig.add_subplot(gs[0 , 0])
            ax.set_title("Figur 1: Identitet per genotyp")

            x_cat = "pident"
            y_cat = "scomname"
            hue_cat = "contig"

            sns.pointplot(
                data=df,
                x=x_cat,
                y=y_cat,
                hue=hue_cat,
                dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
                markers="d", markersize=6, linestyle="none", zorder=10, order = pident_order, legend =False
            )

            sns.stripplot(
                data=df,
                x=x_cat,
                y=y_cat,
                hue=hue_cat,
                dodge=True, alpha=.4, legend=False, jitter = 0.3, order = pident_order
            )

            if n_scomname > 1:
                for i in range(n_scomname):
                    if i % 2 == 1:
                        ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)

            ax.set_ylim(-0.5, n_scomname-0.3)
            ax.grid(True, axis="x", linestyle="--")
            ax.set(xlabel='BLAST identitet (%)', ylabel='Genotyp')

            # Plot 2: bit score
            ax = fig.add_subplot(gs[0 , 1])
            ax.set_title("Figur 2: Bit score per genotyp")

            x_cat = "bitscore"
            y_cat = "scomname"
            hue_cat = "contig"

            sns.pointplot(
                data=df,
                x=x_cat,
                y=y_cat,
                hue=hue_cat,
                dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
                markers="d", markersize=6, linestyle="none", zorder=10, order = bitscore_order, legend = legend_status
            )

            sns.stripplot(
                data=df,
                x=x_cat,
                y=y_cat,
                hue=hue_cat,
                dodge=True, alpha=.4, legend=False, jitter = 0.3, order = bitscore_order
            )
            if n_scomname > 1:
                for i in range(n_scomname):
                    if i % 2 == 1:
                        ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)

            ax.set_ylim(-0.5, n_scomname-0.3)
            ax.grid(True, axis="x", linestyle="--")

            if legend_status == True:
                ax.legend().set_title("Contig")
                sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
            ax.set(xlabel='BLAST bit score', ylabel=' ')

            # Plot 3: coverage
            ax = fig.add_subplot(gs[1, 0])
            ax.set_title("Figur 3: Täckning per contig")

            x_cat = "coverage"
            y_cat = "contig"
            hue_cat = y_cat

            sns.pointplot(
                data=u_contigs,
                x=x_cat,
                y=y_cat,
                hue=hue_cat,
                dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
                markers="d", markersize=6, linestyle="none", order = coverage_order
            )

            if n_contigs > 1:
                for i in range(n_contigs):
                    if i % 2 == 0:
                        ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)

            ax.set_ylim(-0.5, n_contigs-0.3)
            ax.grid(True, axis="x", linestyle="--")
            ax.set(xlabel='Täckning (x)', ylabel='Contig')
            ax.set_xscale('log')

            # Safe axis limits for coverage
            max_coverage = df['coverage'].max()
            if max_coverage > 0:
                ax.set_xlim(left=max(0.1, df['coverage'].min()/2), right=10*max_coverage)
            else:
                ax.set_xlim(left=0.1, right=100)

            # Plot 4: length
            ax = fig.add_subplot(gs[1, 1])
            ax.set_title("Figur 4: Längd per contig")

            x_cat = "qlen"
            y_cat = "contig"
            hue_cat = y_cat

            sns.pointplot(
                data=l_contigs,
                x=x_cat,
                y=y_cat,
                hue=hue_cat,
                dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
                markers="d", markersize=6, linestyle="none", order = length_order, legend = False
            )

            if n_contigs > 1:
                for i in range(n_contigs):
                    if i % 2 == 0:
                        ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)

            ax.set_ylim(-0.5, n_contigs-0.3)
            ax.grid(True, axis="x", linestyle="--")
            ax.set(xlabel='Längd (bp)', ylabel=' ')

            # Safe axis limits for length
            max_length = df['qlen'].max()
            if max_length > 0:
                ax.set_xlim(left=max(1, df['qlen'].min()/2), right=1.1*max_length)
            else:
                ax.set_xlim(left=1, right=1000)

            fig.subplots_adjust(hspace = 20, wspace=20)
            fig.tight_layout()

            # Convert plot to base64 with high resolution
            img_buffer = BytesIO()
            fig.savefig(img_buffer, format='png', dpi=plot_dpi, bbox_inches='tight',
                        facecolor='white', edgecolor='none')
            img_buffer.seek(0)
            img_base64 = base64.b64encode(img_buffer.read()).decode('utf-8')
            img_data_uri = f"data:image/png;base64,{img_base64}"

        except Exception as e:
            print(f"Error generating plots for {file_name}: {e}")
            # Use placeholder image if plot generation fails
            img_data_uri = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=="

        # Read fasta file (use original unfiltered file) - this runs regardless of plot success/failure
        try:
            seq_file_original = f"{output_base}/{ticket_nr}/ev_contig/{file_name.replace('.blast', '.fasta')}"

            if os.path.exists(seq_file_original):
                with open(seq_file_original, 'r') as contigs_file:
                    contigs_content = contigs_file.read()

                if contigs_content:
                    # Extract contig headers for summary
                    contig_headers = extract_contig_headers(contigs_content)
                    contig_count = len(contig_headers)

                    # Create summary (just headers)
                    contigs_summary = '<br>'.join(contig_headers)

                    # Convert full content newlines to HTML breaks for display
                    contigs_content = contigs_content.replace('\n', '<br>')
            else:
                print(f"Warning: FASTA file not found: {seq_file_original}")

        except Exception as e:
            print(f"Error reading FASTA file for {file_name}: {e}")

        # Always return the 4 expected values
        return img_data_uri, contigs_content, contigs_summary, contig_count

    def generate_report(file_path, is_error=False):
        """Generate report using unified template"""
        file_name = os.path.basename(file_path)
        seq_name = file_name.replace('.blast', '')

        if is_error:
            # Use error logic (unfiltered data)
            img_data_uri, contigs_content, contigs_summary, contig_count = generate_error_report_data(file_path)
            varningstext = ''  # No low identity warning for error reports
        else:
            # Use normal logic (filtered data)
            img_data_uri, contigs_content, contigs_summary, contig_count, varningstext = generate_normal_report_data(file_path)

        # Prepare template data
        template_data = {
            'seq_name': seq_name,
            'ticket_nr': ticket_nr,
            'time_stamp': time_stamp,
            'is_error_report': is_error,
            'varningstext': varningstext,
            'img_data_uri': img_data_uri,
            'contigs_content': contigs_content,
            'contigs_summary': contigs_summary,
            'contig_count': contig_count
        }

        # Render unified template
        template = jinja_env.from_string(UNIFIED_TEMPLATE)
        html_content = template.render(**template_data)

        # Create report folder if it doesn't exist
        report_folder = f"{output_base}/{ticket_nr}/report"
        os.makedirs(report_folder, exist_ok=True)

        # Save the HTML content to a file
        with open(f"{report_folder}/{seq_name}.html", "w", encoding='utf-8') as html_file:
            html_file.write(html_content)

    # Main processing - process only the single specified file
    print("\n" + "="*50)
    print(f"Processing file: {blast_file}")
    print("="*50)

    time_stamp = datetime.datetime.now().strftime("%Y-%m-%d")

    try:
        # Try normal report first
        generate_report(blast_file, is_error=False)
        print(f"✅ Success: {os.path.basename(blast_file)}")
        plt.close()
    except (ValueError, FileNotFoundError) as e:
        print(f"⚠️  Normal report failed for {os.path.basename(blast_file)}: {e}")
        print(f"    Generating warning report for {os.path.basename(blast_file)} (using unfiltered contigs)")
        try:
            generate_report(blast_file, is_error=True)
            print(f"⚠️  Warning report generated for {os.path.basename(blast_file)}")
            plt.close()
        except Exception as e2:
            print(f"❌ Failed to generate warning report for {os.path.basename(blast_file)}: {e2}")
            plt.close()
    except Exception as e:
        print(f"❌ Unexpected error for {os.path.basename(blast_file)}: {e}")
        print(f"    Generating warning report for {os.path.basename(blast_file)}")
        try:
            generate_report(blast_file, is_error=True)
            print(f"⚠️  Warning report generated for {os.path.basename(blast_file)}")
            plt.close()
        except Exception as e2:
            print(f"❌ Failed to generate error report for {os.path.basename(blast_file)}: {e2}")
            plt.close()

    # Summary
    print("\n" + "="*50)
    print("PROCESSING COMPLETE")
    print("="*50)
    print(f"📁 Report saved to: {output_base}/{ticket_nr}/report/")
    print("="*50)

if __name__ == "__main__":
    main()