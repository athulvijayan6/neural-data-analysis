SRCDIR = ./src
INCDIR = ./include
BINDIR = ./bin
RLCSINC = ../RLCS/include
RLCSLIB = ../RLCS/lib

SHELL  =  /bin/bash
CC     =  gcc
CFLAGS  =   -Wall

#all: $(BINDIR)/motif_search $(BINDIR)/motif_search_without_group_2p $(BINDIR)/read_pitch_file_sample $(BINDIR)/get_regions_of_interest  $(BINDIR)/sort_spot_result $(BINDIR)/parse_motif_results $(BINDIR)/parse_motif_results_without_group_2p $(BINDIR)/motif_search_only_2p $(BINDIR)/motif_search_only_1p $(BINDIR)/parse_motif_results_only_1p_or_2p  $(BINDIR)/parse_motif_results_across_oneliners $(BINDIR)/motif_search_across_oneliners $(BINDIR)/motif_search_query_region $(BINDIR)/parse_motif_results_query_region $(BINDIR)/parse_motif_filtering $(BINDIR)/join_spot_files $(BINDIR)/gen_ch_wave_file_from_spot_file $(BINDIR)/update_score_with_slope $(BINDIR)/parse_motif_filtering2 $(BINDIR)/motif_search_only_2p_with_filtering $(BINDIR)/motif_search_query_region_with_filtering $(BINDIR)/parse_motif_raga_ID


all: $(BINDIR)/read_pitch_file_sample $(BINDIR)/get_regions_of_interest  $(BINDIR)/sort_spot_result $(BINDIR)/parse_motif_results_top_n_for_each_qmotif $(BINDIR)/parse_motif_results_top_n_for_each_ref $(BINDIR)/parse_motif_results_without_group_2p $(BINDIR)/motif_search_only_2p $(BINDIR)/parse_motif_results_only_1p_or_2p $(BINDIR)/motif_search_across_oneliners $(BINDIR)/motif_search_across_oneliners_with_filtering $(BINDIR)/parse_motif_results_across_oneliners $(BINDIR)/motif_search_query_region $(BINDIR)/parse_motif_results_query_region $(BINDIR)/parse_motif_filtering $(BINDIR)/join_spot_files $(BINDIR)/gen_ch_wave_file_from_spot_file $(BINDIR)/update_score_with_slope $(BINDIR)/parse_motif_filtering2 $(BINDIR)/motif_search_only_2p_with_filtering $(BINDIR)/motif_search_query_region_with_filtering $(BINDIR)/parse_motif_raga_ID $(BINDIR)/yamlwrappertest $(BINDIR)/raga_verification_oneliners $(BINDIR)/motif_search_query_region_2p $(BINDIR)/raga_verification_oneliners_slcs $(BINDIR)/raga_verification_oneliners_rlcs $(BINDIR)/for_plot_slcs $(BINDIR)/raga_verification_oneliners_rlcs_for_znorm $(BINDIR)/raga_verification_oneliners_slcs_for_znorm  $(BINDIR)/raga_verification_oneliners_dtw $(BINDIR)/raga_verification_oneliners_dtw_for_znorm  $(BINDIR)/raga_verification_oneliners_slcs_hard $(BINDIR)/raga_verification_oneliners_slcs_hard_for_znorm  


NAME1=parse_motif_raga_ID
$(BINDIR)/$(NAME1): $(SRCDIR)/$(NAME1).c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/$(NAME1) $(SRCDIR)/$(NAME1).c  -lrlcs -lm

NAME2=update_score_with_slope
$(BINDIR)/$(NAME2): $(SRCDIR)/$(NAME2).c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/$(NAME2) $(SRCDIR)/$(NAME2).c  -lrlcs -lm

NAME3=motif_search_only_2p_with_filtering
$(BINDIR)/$(NAME3): $(SRCDIR)/$(NAME3).c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/$(NAME3) $(SRCDIR)/$(NAME3).c  -lrlcs -lm

NAME4=motif_search_query_region_with_filtering
$(BINDIR)/$(NAME4): $(SRCDIR)/$(NAME4).c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/$(NAME4) $(SRCDIR)/$(NAME4).c  -lrlcs -lm

NAME5=motif_search_across_oneliners_with_filtering
$(BINDIR)/$(NAME5): $(SRCDIR)/$(NAME5).c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/$(NAME5) $(SRCDIR)/$(NAME5).c  -lrlcs -lm

$(BINDIR)/gen_ch_wave_file_from_spot_file: $(SRCDIR)/gen_ch_wave_file_from_spot_file.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/gen_ch_wave_file_from_spot_file $(SRCDIR)/gen_ch_wave_file_from_spot_file.c  -lrlcs -lm

$(BINDIR)/join_spot_files: $(SRCDIR)/join_spot_files.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/join_spot_files $(SRCDIR)/join_spot_files.c  -lrlcs -lm

$(BINDIR)/motif_search_query_region: $(SRCDIR)/motif_search_query_region.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/motif_search_query_region $(SRCDIR)/motif_search_query_region.c  -lrlcs -lm

$(BINDIR)/motif_search: $(SRCDIR)/motif_search.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/motif_search $(SRCDIR)/motif_search.c  -lrlcs -lm

$(BINDIR)/motif_search_across_oneliners: $(SRCDIR)/motif_search_across_oneliners.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/motif_search_across_oneliners $(SRCDIR)/motif_search_across_oneliners.c  -lrlcs -lm

$(BINDIR)/motif_search_only_1p: $(SRCDIR)/motif_search_only_1p.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/motif_search_only_1p $(SRCDIR)/motif_search_only_1p.c  -lrlcs -lm

$(BINDIR)/parse_motif_filtering2: $(SRCDIR)/parse_motif_filtering2.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/parse_motif_filtering2 $(SRCDIR)/parse_motif_filtering2.c -lrlcs -lm

$(BINDIR)/parse_motif_filtering: $(SRCDIR)/parse_motif_filtering.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/parse_motif_filtering $(SRCDIR)/parse_motif_filtering.c -lrlcs -lm

$(BINDIR)/parse_motif_results_query_region: $(SRCDIR)/parse_motif_results_query_region.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/parse_motif_results_query_region $(SRCDIR)/parse_motif_results_query_region.c -lrlcs -lm

$(BINDIR)/parse_motif_results_only_1p_or_2p: $(SRCDIR)/parse_motif_results_only_1p_or_2p.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/parse_motif_results_only_1p_or_2p $(SRCDIR)/parse_motif_results_only_1p_or_2p.c  -lrlcs -lm

$(BINDIR)/parse_motif_results_across_oneliners: $(SRCDIR)/parse_motif_results_across_oneliners.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/parse_motif_results_across_oneliners $(SRCDIR)/parse_motif_results_across_oneliners.c  -lrlcs -lm

$(BINDIR)/motif_search_only_2p: $(SRCDIR)/motif_search_only_2p.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/motif_search_only_2p $(SRCDIR)/motif_search_only_2p.c  -lrlcs -lm

$(BINDIR)/motif_search_without_group_2p: $(SRCDIR)/motif_search_without_group_2p.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/motif_search_without_group_2p $(SRCDIR)/motif_search_without_group_2p.c  -lrlcs -lm

$(BINDIR)/read_pitch_file_sample: $(SRCDIR)/read_pitch_file_sample.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/read_pitch_file_sample $(SRCDIR)/read_pitch_file_sample.c  -lrlcs -lm

$(BINDIR)/get_regions_of_interest: $(SRCDIR)/get_regions_of_interest.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/get_regions_of_interest $(SRCDIR)/get_regions_of_interest.c  -lrlcs -lm

$(BINDIR)/sort_spot_result: $(SRCDIR)/sort_spot_result.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/sort_spot_result $(SRCDIR)/sort_spot_result.c -lrlcs -lm


$(BINDIR)/parse_motif_results_top_n_for_each_qmotif: $(SRCDIR)/parse_motif_results_top_n_for_each_qmotif.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/parse_motif_results_top_n_for_each_qmotif $(SRCDIR)/parse_motif_results_top_n_for_each_qmotif.c  -lrlcs -lm

$(BINDIR)/parse_motif_results_top_n_for_each_ref: $(SRCDIR)/parse_motif_results_top_n_for_each_ref.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/parse_motif_results_top_n_for_each_ref $(SRCDIR)/parse_motif_results_top_n_for_each_ref.c  -lrlcs -lm

$(BINDIR)/parse_motif_results_without_group_2p: $(SRCDIR)/parse_motif_results_without_group_2p.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/parse_motif_results_without_group_2p $(SRCDIR)/parse_motif_results_without_group_2p.c  -lrlcs -lm

$(BINDIR)/yamlwrappertest: $(SRCDIR)/yamlwrappertest.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/yamlwrappertest $(SRCDIR)/yamlwrappertest.c  -lrlcs -lm


$(BINDIR)/raga_verification_oneliners: $(SRCDIR)/raga_verification_oneliners.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/raga_verification_oneliners $(SRCDIR)/raga_verification_oneliners.c  -lrlcs -lm


$(BINDIR)/motif_search_query_region_2p: $(SRCDIR)/motif_search_query_region_2p.c
	$(CC) $(CFFLAGS) -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/motif_search_query_region_2p $(SRCDIR)/motif_search_query_region_2p.c  -lrlcs -lm

$(BINDIR)/raga_verification_oneliners_slcs: $(SRCDIR)/raga_verification_oneliners_slcs.c
	$(CC) $(CFFLAGS) -fopenmp -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/raga_verification_oneliners_slcs $(SRCDIR)/raga_verification_oneliners_slcs.c  -lrlcs -lm


$(BINDIR)/raga_verification_oneliners_rlcs: $(SRCDIR)/raga_verification_oneliners_rlcs.c
	$(CC) $(CFFLAGS) -fopenmp -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/raga_verification_oneliners_rlcs $(SRCDIR)/raga_verification_oneliners_rlcs.c  -lrlcs -lm

$(BINDIR)/raga_verification_oneliners_rlcs_for_znorm: $(SRCDIR)/raga_verification_oneliners_rlcs_for_znorm.c
	$(CC) $(CFFLAGS) -fopenmp -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/raga_verification_oneliners_rlcs_for_znorm $(SRCDIR)/raga_verification_oneliners_rlcs_for_znorm.c  -lrlcs -lm

$(BINDIR)/raga_verification_oneliners_slcs_for_znorm: $(SRCDIR)/raga_verification_oneliners_slcs_for_znorm.c
	$(CC) $(CFFLAGS) -fopenmp -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/raga_verification_oneliners_slcs_for_znorm $(SRCDIR)/raga_verification_oneliners_slcs_for_znorm.c  -lrlcs -lm


$(BINDIR)/raga_verification_oneliners_dtw: $(SRCDIR)/raga_verification_oneliners_dtw.c
	$(CC) $(CFFLAGS) -fopenmp -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/raga_verification_oneliners_dtw $(SRCDIR)/raga_verification_oneliners_dtw.c  -lrlcs -lm


$(BINDIR)/raga_verification_oneliners_dtw_for_znorm: $(SRCDIR)/raga_verification_oneliners_dtw_for_znorm.c
	$(CC) $(CFFLAGS) -fopenmp -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/raga_verification_oneliners_dtw_for_znorm $(SRCDIR)/raga_verification_oneliners_dtw_for_znorm.c  -lrlcs -lm

$(BINDIR)/raga_verification_oneliners_slcs_hard: $(SRCDIR)/raga_verification_oneliners_slcs_hard.c
	$(CC) $(CFFLAGS) -fopenmp -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/raga_verification_oneliners_slcs_hard $(SRCDIR)/raga_verification_oneliners_slcs_hard.c  -lrlcs -lm


$(BINDIR)/raga_verification_oneliners_slcs_hard_for_znorm: $(SRCDIR)/raga_verification_oneliners_slcs_hard_for_znorm.c
	$(CC) $(CFFLAGS) -fopenmp -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/raga_verification_oneliners_slcs_hard_for_znorm $(SRCDIR)/raga_verification_oneliners_slcs_hard_for_znorm.c  -lrlcs -lm



$(BINDIR)/for_plot_slcs: $(SRCDIR)/for_plot_slcs.c
	$(CC) $(CFFLAGS) -fopenmp -L $(RLCSLIB) -I $(RLCSINC) -o $(BINDIR)/for_plot_slcs $(SRCDIR)/for_plot_slcs.c  -lrlcs -lm



clean:
	rm -rf  $(BINDIR)/*
