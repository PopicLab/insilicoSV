from collections import defaultdict
from contextlib import closing
import subprocess
from pysam import VariantFile

sr_path = '/data/projects/insilico/igv_figures/reads/sim_sr.dwgsim.minimap2.sorted.bam'
lr_path = '/data/projects/insilico/igv_figures/reads/sim_lr.minimap2.sorted.bam'
chr = 'chr1'
output_path ='/data/projects/insilico/igv_figures/figures/'
vcf_path = '/data/projects/insilico/igv_figures/output/sim.vcf'

def call_samplot(filename):
    recs = defaultdict(list)
    num_simple_sv = 0
    with closing(VariantFile(filename)) as vcf:
        bnd = []
        for vcf_rec in vcf.fetch():
            vcf_info = dict(vcf_rec.info)
            type = vcf_info['SVTYPE']
            positions = [vcf_rec.start, vcf_rec.stop]
            letter = vcf_info['SOURCE_LETTER']
            target = None
            if positions in bnd: continue
            bnd.append(positions)
            if 'TARGET' in vcf_info:
                target = vcf_info['TARGET']
            if 'PARENT_SVID' in vcf_info:
                recs[vcf_info['PARENT_SVID']].append([positions, target, vcf_info['PARENT_SVTYPE']+'_'+vcf_info['PARENT_SVID'], letter])
            else:
                recs[str(num_simple_sv)].append([positions, target, type + '_' + vcf_rec.id, letter])
                num_simple_sv += 1
    for parent_id, sv_recs in recs.items():
        sv_breakends = []
        letters = {}
        targets = []
        for sv_rec in sv_recs:
            breakends, target, sv_type, letter = sv_rec
            letters[letter] = breakends[:1]
            sv_breakends += breakends
            if target is not None:
                targets.append(target)
        sv_breakends = sorted(sv_breakends)
        output_file = output_path + sv_type + '.png'
        targets = [target for target in targets if target not in sv_breakends]
        starts = [min([str(breakend) for idx, breakend in enumerate(sv_breakends) if idx % 2 == 0] + [str(target) for target in targets])]
        ends = [max([str(breakend) for idx, breakend in enumerate(sv_breakends) if idx % 2 == 1] + [str(target +1) for target in targets])]

        subprocess.run(['samplot', 'plot', '-n', 'Short-reads Long-reads', '-b', sr_path, lr_path, '-s'] + starts +
                        ['-e'] + ends + ['-c', chr, '-t', sv_type, '--separate_mqual', '1', '--include_mqual', '0',
                        '--legend_fontsize', '10', '-o', output_file])
call_samplot(vcf_path)