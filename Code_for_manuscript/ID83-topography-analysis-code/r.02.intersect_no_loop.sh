#!/bin/bash
nohup intersectBed -a /public/data/Topography_analysis_2026/results/ID83/01_createBed/01.rtGroup.bed -b /public/data/Topography_analysis_2026/results0322/ID83/01_createBed/03.signature.bed -wa -wb > /public/data/Topography_analysis_2026/results0322/ID83/02_intersectBed/01.RelicateTime.TransStrand.DNA_Region.intersect_res.txt &
nohup intersectBed -a /public/data/Topography_analysis/reference/ENCODE/repli_strand_info_new.bed -b /public/data/Topography_analysis_2026/results0322/ID83/01_createBed/03.signature.bed -wa -wb > /public/data/Topography_analysis_2026/results0322/ID83/02_intersectBed/02.RepliStrand.intersect_res.txt &

nohup intersectBed -a /public/data/Topography_analysis_2026/results/ID83/01_createBed/01.rtGroup.bed -b /public/data/Topography_analysis_2026/results0322/ID83/01_createBed/04.SimulatedData.bed -wa -wb > /public/data/Topography_analysis_2026/results0322/ID83/02_intersectBed/03.RelicateTime.simulated.intersect_res.txt &
nohup intersectBed -a /public/data/Topography_analysis/reference/ENCODE/repli_strand_info_new.bed -b /public/data/Topography_analysis_2026/results0322/ID83/01_createBed/04.SimulatedData.bed -wa -wb > /public/data/Topography_analysis_2026/results0322/ID83/02_intersectBed/04.RepliStrand.simulated.intersect_res.txt &

