import pandas as pd
import numpy as np
import os
import pickle
import ast

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import numpy as np
import math

import warnings

from matplotlib.backends.backend_pdf import PdfPages

os.chdir('/data3/projects/2025_Antibiotics/YSH/')

### 2 ###

def make_input(datadir, pids):
    sev_dict = {}
    sev_dict_filtered = {}

    for pid in pids:
        try:
            # 1. SeverityScore 불러오기
            sev = pd.read_csv(f"{datadir}/{pid}/SeverityScore.csv")

            # 2. Laboratory 데이터 불러오기 및 병합
            lab = pd.read_csv(f"{datadir}/{pid}/Laboratory_processed.csv")
            sev = pd.merge(sev, lab, on='Date', how='left')

            # 3. Medication 불러오기
            med = pd.read_csv(f"{datadir}/{pid}/Medication.csv")
            med_filtered = med[med['Date'].isin(sev['Date'])]
            med_filtered = med_filtered.loc[:, ~med_filtered.columns.str.endswith('_dose')]

            # 4. Medication 관련 열 추가
            sev['med_cnt'] = 0
            sev['med_list'] = ""

            # 5. 날짜별로 약물 정보 병합
            for _, row in med_filtered.iterrows():
                cur_date = row['Date']
                cur_meds_raw = row.iloc[1:]
                cur_meds_clean = cur_meds_raw.dropna().tolist()

                med_freq = {}
                cur_meds = []
                for med in cur_meds_clean:
                    if med not in med_freq:
                        med_freq[med] = 1
                        cur_meds.append(med)
                    else:
                        med_freq[med] += 1
                        cur_meds.append(f"{med}_{med_freq[med]}")

                cur_med_cnt = len(cur_meds)
                cur_med_str = ";".join(cur_meds)

                sev_idx = sev[sev['Date'] == cur_date].index
                if len(sev_idx) > 0:
                    sev.loc[sev_idx, 'med_cnt'] = cur_med_cnt
                    sev.loc[sev_idx, 'med_list'] = cur_med_str

            zero_cnt = (sev['med_cnt'] == 0).sum()

        except FileNotFoundError:
            sev['med_cnt'] = ""
            sev['med_list'] = ""
            zero_cnt = "N/A (no med file)"

        if sev.shape[0] == zero_cnt:
            print(pid)

        if zero_cnt == "N/A (no med file)":
            print(pid, zero_cnt)

        sev_dict[pid] = sev

        # Remove med cnt 0
        if 'med_cnt' not in sev.columns:
            print(f"{pid} - Error: 'med_cnt' column doesn't exist")
            continue
        sev['med_cnt'] = pd.to_numeric(sev['med_cnt'], errors='coerce')

        if sev['med_cnt'].isna().all():
            continue

        if sev['med_cnt'].fillna(0).eq(0).all():
            continue

        sev_dict_filtered[pid] = sev
        
    return sev_dict_filtered
    
def add_strain(sev_dict_filtered):
    strain_path = '/data3/projects/2025_Antibiotics/data/병원체자원은행 균주현황(2014-2024.06)_Sepsis.xlsx'
    strain_df = pd.read_excel(strain_path)

    valid_columns = strain_df.iloc[0].dropna().index
    strain_legend = strain_df.loc[:1, valid_columns].copy()

    strain_df.columns = strain_df.iloc[1]
    strain_df = strain_df.drop(index=[0, 1])
    strain_df = strain_df[['접수일', '등록번호', '균']]
    strain_df['등록번호'] = strain_df['등록번호'].astype(int)

    valid_ids = [int(k) for k in sev_dict_filtered.keys()]
    #strain_df['접수번호'] = strain_df['접수번호'].astype(str).astype(int)
    strain_df = strain_df[strain_df['등록번호'].isin(valid_ids)].reset_index(drop=True)

    strain_df = strain_df.drop_duplicates(subset=['균', '등록번호'], keep='first').reset_index(drop=True)
    strain_df['접수일'] = pd.to_datetime(strain_df['접수일'])
    strain_df['접수일'] = strain_df['접수일'].dt.strftime('%Y-%m-%d')

    # 등록번호 컬럼을 str로 변환 (딕셔너리 key와 비교하기 위함)
    strain_df['등록번호'] = strain_df['등록번호'].astype(str)

    pid_without_strains = []   

    for cur_key, cur_sev_df in sev_dict_filtered.items():
        cur_key_str = str(cur_key)

        # 등록번호가 현재 key인 행들 필터링
        cur_df = strain_df[strain_df['등록번호'] == cur_key_str]

        if cur_df.empty:
            print(cur_key)
            continue

        # 날짜 변환
        cur_sev_df['Date'] = pd.to_datetime(cur_sev_df['Date'])
        cur_df['접수일'] = pd.to_datetime(cur_df['접수일'])

        for _, row in cur_df.iterrows():
            cur_date = row['접수일']
            cur_strain = row['균']

            # 접수일과 같거나 이후인 첫 행의 인덱스 찾기
            matched_idx = cur_sev_df[cur_sev_df['Date'] >= cur_date].index.min()

            if pd.isna(matched_idx):
                continue  # 매칭된 날짜가 없다면 skip

            # strain 컬럼이 없다면 빈 리스트 생성
            if 'strain' not in cur_sev_df.columns:
                cur_sev_df['strain'] = [[] for _ in range(len(cur_sev_df))]

            # strain 컬럼이 리스트 형식이 아니면 변환
            if not isinstance(cur_sev_df.at[matched_idx, 'strain'], list):
                cur_sev_df.at[matched_idx, 'strain'] = []

            # 중복 방지를 원할 경우 다음 줄에 조건 추가 가능
            cur_sev_df.at[matched_idx, 'strain'].append(cur_strain)

        # strain 컬럼이 없는 경우 스킵
        if 'strain' not in cur_sev_df.columns:
            pid_without_strains.append(cur_key_str)
            continue

        cur_sev_df = cur_sev_df.reset_index(drop=True)

        last_strain = []  # 최근에 발견된 strain 리스트
        for i in range(len(cur_sev_df)):
            current_strain = cur_sev_df.at[i, 'strain']

            if isinstance(current_strain, list) and len(current_strain) > 0:
                # 비어있지 않은 strain 리스트 발견 → 이를 저장
                last_strain = current_strain
            elif isinstance(current_strain, list) and len(current_strain) == 0:
                # 비어있다면 → 가장 최근의 strain 리스트를 할당
                cur_sev_df.at[i, 'strain'] = last_strain.copy()

        # 각 strain 리스트에서 중복 제거
        cur_sev_df['strain'] = cur_sev_df['strain'].apply(
            lambda x: list(set(x)) if isinstance(x, list) else x
        )

        # 다시 딕셔너리에 반영
        sev_dict_filtered[cur_key] = cur_sev_df
        
    return sev_dict_filtered, pid_without_strains

    
### 3 ###

def make_sev_dict(med, indir, outdir):
    with open(f"{indir}/Input.pkl", 'rb') as f:
        sev_dict_filtered = pickle.load(f)
    
    if med != "total":
        for pid, cur_df in sev_dict_filtered.items():
            try:
                cur_df['cur_med_list'] = cur_df['med_list'].fillna('').apply(lambda x: x.split(';'))
                cur_df['med_cnt'] = cur_df['cur_med_list'].apply(
                    lambda meds: sum([1 for item in meds if item.startswith(med)])
                )
                #cur_df = cur_df[['Date', 'NEWS', 'med_cnt', 'strain']]
                #print(cur_df.shape)
                sev_dict_filtered[pid] = cur_df
            except KeyError as e:
                print(f"[{pid}] 건너뜀 - 누락 컬럼: {e}")
                continue
    else:
        for pid, cur_df in sev_dict_filtered.items():
            try:
                #cur_df = cur_df[['Date', 'NEWS', 'med_cnt', 'strain']]
                sev_dict_filtered[pid] = cur_df
            except KeyError as e:
                print(f"[{pid}] 건너뜀 - 누락 컬럼: {e}")
                continue

    with open(f"{outdir}/{med}.pkl", 'wb') as f:
        pickle.dump(sev_dict_filtered, f)
    
    return sev_dict_filtered

def summarize_med_cnt(med_cnt_series):
    info_list = []
    if med_cnt_series.empty:
        return info_list

    cur_type = 'm' if med_cnt_series.iloc[0] > 0 else 'n'
    count = 1

    for prev, curr in zip(med_cnt_series[:-1], med_cnt_series[1:]):
        curr_type = 'm' if curr > 0 else 'n'
        if curr_type == cur_type:
            count += 1
        else:
            info_list.append(f"{count}{cur_type}")
            cur_type = curr_type
            count = 1

    info_list.append(f"{count}{cur_type}")
    return info_list

def make_timecourse(indir, outdir, med):
    with open(f"{indir}/{med}.pkl", 'rb') as f:
        sev_dict_filtered = pickle.load(f)

    timecourse_list = []

    for pid, sev in sev_dict_filtered.items():
        if 'med_cnt' not in sev.columns or sev['med_cnt'].isnull().all():
            print(f"{pid} - Error: 'med_cnt' column doesn't exist")
            continue

        sev['med_cnt'] = pd.to_numeric(sev['med_cnt'], errors='coerce').fillna(0).astype(int)

        info_list = summarize_med_cnt(sev['med_cnt'])

        timecourse_list.append({
            'pid': pid,
            'days': len(sev),
            'info': info_list
        })

    timecourse = pd.DataFrame(timecourse_list)
    timecourse.to_csv(f"{outdir}/{med}.csv", index=False)
    
    return timecourse

def filter_info(info_list):
    has_m_prefix = any(item.startswith('m') for item in info_list)
    has_m_suffix = any(item.endswith('m') for item in info_list)
    has_n_suffix = any(item.endswith('n') for item in info_list)

    if has_m_prefix:
        return False
    if all(item.endswith('m') for item in info_list):
        return False
    if all(item.endswith('n') for item in info_list):
        return False
    return True

def convert_info_with_repeats(info_list):
    result = ''
    for item in info_list:
        if item.endswith('n'):
            num = int(item[:-1])
            result += '_' * num
        elif item.endswith('m'):
            num = int(item[:-1])
            result += 'm' * num
    return result

def find_valid_start_end(s):
    valid_starts = []
    valid_ends = []
    for i in range(len(s) - 3):
        if s[i:i+4] == '___m':
            end = i + 9
            if end < len(s):
                valid_starts.append(i)
                valid_ends.append(end)
    return pd.Series({'start_idx': valid_starts, 'end_idx': valid_ends})

#def find_valid_start_end(s):
#    valid_starts = []
#    valid_ends = []
#    for i in range(len(s) - 7): 
#        if s[i:i+8] == '_______m':
#            end = i + 13  
#            if end < len(s):
#                valid_starts.append(i)
#                valid_ends.append(end)
#    return pd.Series({'start_idx': valid_starts, 'end_idx': valid_ends})

def make_sev_idx(indir, outdir, med):
    timecourse = pd.read_csv(f"{indir}/{med}.csv")

    if isinstance(timecourse['info'].iloc[0], str):
        timecourse['info'] = timecourse['info'].apply(ast.literal_eval)

    timecourse = timecourse[timecourse['info'].apply(filter_info)].reset_index(drop=True)
    timecourse['info_timecourse'] = timecourse['info'].apply(convert_info_with_repeats)
    
    if not timecourse.empty:
        timecourse[['start_idx', 'end_idx']] = timecourse['info_timecourse'].apply(find_valid_start_end)
        timecourse = timecourse[~((timecourse['start_idx'].apply(len) == 0) & (timecourse['end_idx'].apply(len) == 0))].reset_index(drop=True)
    
    timecourse.to_csv(f"{outdir}/{med}.csv", index=False)

    return timecourse


def make_res_dict(indir, indexdir, outdir, med):
    
    with open(f"{indir}/{med}.pkl", 'rb') as f:
        sev_dict_filtered = pickle.load(f)

    sev_idx = pd.read_csv(f"{indexdir}/{med}.csv")
    
    res_dict = {}

    if not sev_idx.empty:  
        if isinstance(sev_idx['start_idx'].iloc[0], str):
            sev_idx['start_idx'] = sev_idx['start_idx'].apply(ast.literal_eval)

        if isinstance(sev_idx['end_idx'].iloc[0], str):
            sev_idx['end_idx'] = sev_idx['end_idx'].apply(ast.literal_eval)

        for _, row in sev_idx.iterrows():
            cur_pid = row['pid']
            cur_start_list = row['start_idx']
            cur_end_list = row['end_idx']
            cur_df = sev_dict_filtered[str(cur_pid)]

            for i in range(len(cur_start_list)):
                cur_start_idx = cur_start_list[i]
                cur_end_idx = cur_end_list[i]
                cur_filtered_df = cur_df.iloc[cur_start_idx:cur_end_idx+1].copy()
                key = f"{cur_pid}_{i}"
                res_dict[key] = cur_filtered_df

    with open(f"{outdir}/{med}.pkl", 'wb') as f:
        pickle.dump(res_dict, f)
        
    return res_dict

def make_sequence(med, indir, outdir):
    sevdir = f'data_{dtype}/temp/sev_dict'
    sev_dict = make_sev_dict(med, indir, sevdir)
    
    timecourse_dir = f'data_{dtype}/temp/timecourse'
    timecourse = make_timecourse(sevdir, timecourse_dir, med)
    
    idx_dir = f'data_{dtype}/temp/sev_idx'
    sev_idx = make_sev_idx(timecourse_dir, idx_dir, med)
    
    res_dict = make_res_dict(sevdir, idx_dir, outdir, med)
    print(f'Sequence saved as {outdir}/{med}.pkl')
    
    
    
    
    
    
### 4 ###


### temp ###

def find_optimal_k_with_plot(X, min_k, max_k, plot=True, allow_weak_peaks=True, score_ratio_threshold=0.95):
    if len(X) < min_k:
        print(f"샘플 수(n={len(X)})가 {min_k}보다 적어 k={len(X)-1}으로 설정")
        return len(X)

    if len(X) == min_k:
        print(f"샘플 수(n={len(X)})가 {min_k}으로 k={len(X)-1}으로 설정")
        return min_k
    
    k_values = list(range(min_k, max_k + 1))
    
    if len(X) < max_k:
        print(f"샘플 수(n={len(X)})가 {max_k}보다 적어 k={len(X)-1}까지 탐색")
        k_values = list(range(min_k, len(X)))
    
    if len(X) == max_k:
        print(f"샘플 수(n={len(X)})가 {max_k}으로 k={len(X)-1}까지 탐색")
        k_values = list(range(min_k, len(X)))
    
    silhouette_scores = []

    for k in k_values:
        kmeans = KMeans(n_clusters=k, random_state=42)
        labels = kmeans.fit_predict(X)
        score = silhouette_score(X, labels, metric="euclidean")
        silhouette_scores.append(score)

    max_score = max(silhouette_scores)
    best_k = k_values[np.argmax(silhouette_scores)]

    score_diffs = np.diff(silhouette_scores)
    peak_candidates = []

    for i in range(1, len(score_diffs)):
        if score_diffs[i] > 0 and score_diffs[i - 1] < 0:
            local_k = k_values[i + 1]
            local_score = silhouette_scores[i + 1]

            if allow_weak_peaks:
                peak_candidates.append((local_k, local_score))
            elif local_score >= score_ratio_threshold * max_score:
                peak_candidates.append((local_k, local_score))

    if peak_candidates:
        optimal_k = max(peak_candidates, key=lambda x: x[1])[0]
    else:
        optimal_k = best_k

    if plot:
        plt.figure(figsize=(10, 6))
        plt.plot(k_values, silhouette_scores, marker='o')
        plt.axvline(optimal_k, color='red', linestyle='--', label=f"Optimal k = {optimal_k}")
        plt.title("Silhouette Scores by Number of Clusters (KMeans)")
        plt.xlabel("Number of clusters (k)")
        plt.ylabel("Silhouette Score (Euclidean)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return optimal_k

def perform_kmeans(indir, sevdir, outdir, med):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        
        with open(f"{indir}/{med}.pkl", 'rb') as f:
            res_dict = pickle.load(f)

        news_dict = {pid: df['NEWS'].tolist() for pid, df in res_dict.items()}
        news_df = pd.DataFrame(news_dict)
        
        if news_df.empty:
            cluster_df = news_df
            print("샘플 수 0개")

        if not news_df.empty:  
            X = news_df.T.values

            optimal_k = find_optimal_k_with_plot(X, min_k=30, max_k=40)

            kmeans = KMeans(n_clusters=optimal_k, random_state=42)
            kmeans.fit(X)
            labels = kmeans.labels_

            cluster_df = news_df.T
            cluster_df = cluster_df.rename(columns=lambda x: f"t{int(x)+1}")
            cluster_df['cluster'] = labels

        cluster_df.to_csv(f"{outdir}/cluster/{med}.csv")

    return cluster_df

def draw_clusterwise_plot(outdir, used_meds):
    
    with PdfPages(f"{outdir}/clusterwise_plot.pdf") as pdf:
        for med in used_meds:

            cluster_df = pd.read_csv(f"{outdir}/cluster/{med}.csv", index_col=0)

            N = cluster_df.shape[1] - 1 if 'cluster' in cluster_df.columns else cluster_df.shape[0]
            time_points = [f"t{i+1}" for i in range(N)]

            clusters = sorted(cluster_df['cluster'].unique())
            n_clusters = len(clusters)

            n_cols = 8
            n_rows = math.ceil(n_clusters / n_cols)

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.5 * n_rows), sharey=True)
            axes = axes.flatten()

            for idx, cluster_id in enumerate(clusters):
                ax = axes[idx]

                cluster_data = cluster_df[cluster_df['cluster'] == cluster_id][time_points]

                norm_rows = []

                for _, row in cluster_data.iterrows():
                    values = row.values.astype(float)

                    norm = np.linalg.norm(values)
                    if norm != 0:
                        norm_values = values / norm
                    else:
                        norm_values = np.zeros_like(values)

                    norm_rows.append(norm_values)
                    ax.plot(time_points, norm_values, color='gray', alpha=0.4)

                norm_array = np.array(norm_rows)
                mean_series = norm_array.mean(axis=0)
                std_series = norm_array.std(axis=0)

                first_half_mean = mean_series[3:4].mean()   # t4
                second_half_mean = mean_series[4:10].mean()  # t5~t10

                if first_half_mean > second_half_mean:
                    mean_color = 'blue'
                else:
                    mean_color = 'blue'

                ax.plot(time_points, mean_series, color=mean_color, linewidth=2, label='Mean')

                ax.fill_between(time_points,
                                mean_series - std_series,
                                mean_series + std_series,
                                color=mean_color,
                                alpha=0.2,
                                label='±1 std')

                if 't8' in time_points:
                    t4_index = time_points.index('t8')
                    ax.axvline(x=time_points[t4_index], color='red', linestyle='--', linewidth=1)

                ax.set_title(f'Cluster {cluster_id} (n={len(cluster_data)})')
                ax.set_xlabel("Time Point")
                ax.set_ylabel("L2 Normalized NEWS")
                ax.set_ylim(0, 1.1)

            plt.suptitle(f"{med} - {cluster_df.shape[0]} samples / k={n_clusters}", fontsize=18)
            plt.tight_layout(rect=[0, 0, 1, 0.96])

            # 현재 figure를 pdf에 저장
            pdf.savefig(fig)
            plt.close(fig)  # 메모리 절약 위해 닫기

def draw_clusterwise_plot_v2(outdir, used_meds):
    with PdfPages(f"{outdir}/clusterwise_plot_v2.pdf") as pdf:
        for med in used_meds:
            cluster_df = pd.read_csv(f"{outdir}/cluster/{med}.csv", index_col=0)

            N = cluster_df.shape[1] - 1 if 'cluster' in cluster_df.columns else cluster_df.shape[0]
            time_points = [f"t{i+1}" for i in range(N)]

            clusters = sorted(cluster_df['cluster'].unique())
            n_clusters = len(clusters)

            n_cols = 8
            n_rows = math.ceil(n_clusters / n_cols)

            fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.5 * n_rows), sharey=True)
            axes = axes.flatten()

            for idx, cluster_id in enumerate(clusters):
                ax = axes[idx]

                cluster_data = cluster_df[cluster_df['cluster'] == cluster_id][time_points]

                raw_rows = []

                for _, row in cluster_data.iterrows():
                    values = row.values.astype(float)

                    raw_rows.append(values)
                    ax.plot(time_points, values, color='gray', alpha=0.4)

                raw_array = np.array(raw_rows)
                mean_series = raw_array.mean(axis=0)
                std_series = raw_array.std(axis=0)

                first_half_mean = mean_series[3:4].mean()   # t4
                second_half_mean = mean_series[4:10].mean()  # t5~t10

                mean_color = 'blue'  # 일단 조건 분기 없어도 됨

                ax.plot(time_points, mean_series, color=mean_color, linewidth=2, label='Mean')

                ax.fill_between(time_points,
                                mean_series - std_series,
                                mean_series + std_series,
                                color=mean_color,
                                alpha=0.2,
                                label='±1 std')

                if 't4' in time_points:
                    t4_index = time_points.index('t4')
                    ax.axvline(x=time_points[t4_index], color='red', linestyle='--', linewidth=1)

                ax.set_title(f'Cluster {cluster_id} (n={len(cluster_data)})')
                ax.set_xlabel("Time Point")
                ax.set_ylabel("Raw NEWS")
                #ax.set_ylim(0, 13)
                ax.set_ylim(bottom=0)

            plt.suptitle(f"{med} - {cluster_df.shape[0]} samples / k={n_clusters}", fontsize=18)
            plt.tight_layout(rect=[0, 0, 1, 0.96])

            pdf.savefig(fig)
            plt.close(fig)
