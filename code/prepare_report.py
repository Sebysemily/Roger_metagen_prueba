#!/usr/bin/env python3
"""
Prepare report directories and copy images.
Usage: python prepare_report.py --html-dir DIR --samples "sample1 sample2" --pre-ybl "file1 file2" ...
"""

import os
import shutil
import argparse
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(description='Prepare QC report directories and images')
    parser.add_argument('--html-dir', required=True, help='HTML output directory')
    parser.add_argument('--samples', required=True, help='Space-separated sample names')
    parser.add_argument('--pre-ybl', required=True, help='Space-separated pre-QC Yield_By_Length.png files')
    parser.add_argument('--pre-lvq', required=True, help='Space-separated pre-QC LengthvsQualityScatterPlot_dot.png files')
    parser.add_argument('--post-ybl', required=True, help='Space-separated post-QC Yield_By_Length.png files')
    parser.add_argument('--post-lvq', required=True, help='Space-separated post-QC LengthvsQualityScatterPlot_dot.png files')
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    # Parse inputs
    html_dir = Path(args.html_dir)
    samples = args.samples.split()
    
    # Parse file lists (they come as space-separated strings)
    pre_ybl_files = args.pre_ybl.split()
    pre_lvq_files = args.pre_lvq.split()
    post_ybl_files = args.post_ybl.split()
    post_lvq_files = args.post_lvq.split()
    
    print(f"📁 Preparing report in: {html_dir}")
    print(f"🧬 Samples ({len(samples)}): {', '.join(samples)}")
    
    # Validate input lengths match
    expected_len = len(samples)
    file_lists = {
        'pre_ybl': pre_ybl_files,
        'pre_lvq': pre_lvq_files,
        'post_ybl': post_ybl_files,
        'post_lvq': post_lvq_files
    }
    
    for name, files in file_lists.items():
        if len(files) != expected_len:
            print(f"⚠️  Warning: {name} has {len(files)} files, expected {expected_len}")
    
    # Create main directory
    html_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each sample
    for i, sample in enumerate(samples):
        sample_img_dir = html_dir / "images" / sample
        sample_img_dir.mkdir(parents=True, exist_ok=True)
        
        print(f"\n📦 Processing sample: {sample}")
        
        # Copy images with error handling
        files_to_copy = [
            (pre_ybl_files[i] if i < len(pre_ybl_files) else None, "Yield_By_Length.png"),
            (pre_lvq_files[i] if i < len(pre_lvq_files) else None, "LengthvsQualityScatterPlot_dot.png"),
            (post_ybl_files[i] if i < len(post_ybl_files) else None, "Yield_By_Length.png"),
            (post_lvq_files[i] if i < len(post_lvq_files) else None, "LengthvsQualityScatterPlot_dot.png")
        ]
        
        for src, dst in files_to_copy:
            if src and os.path.exists(src):
                shutil.copy2(src, sample_img_dir / dst)
                print(f"  ✅ Copied: {dst}")
            else:
                print(f"  ⚠️  Missing: {dst}")
    
    print(f"\n🎉 Preparation complete!")
    print(f"   Report will be at: {html_dir}/01_qc.html")
    print(f"   Images at: {html_dir}/images/")
    
    # Create a simple README in the images directory
    readme = html_dir / "images" / "README.txt"
    readme.write_text(f"""Images organized by sample:
- Each sample has its own subdirectory
- Pre and post-QC images have the same filename
- Images are referenced as: images/{{sample}}/filename.png

Samples: {', '.join(samples)}
Generated: {os.path.basename(__file__)}
""")

if __name__ == "__main__":
    main()
