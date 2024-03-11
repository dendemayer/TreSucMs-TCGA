from PyPDF2 import PdfMerger
import os

"""
taking the list of pdfs which shall be merged and the output path where the
final merged pdf shall be saved to
"""

def merge_pdfs(pdfs_to_merge, out_path):
    slice_length = 100
    mod = len(pdfs_to_merge) % slice_length
    pdfs_to_merge_lists = []
    start_ind = 0
    slices = int(len(pdfs_to_merge) / slice_length)
    if len(pdfs_to_merge) > slice_length:
        # divide the pdfs to merge in slice_lengther steps, get the rest:
        for i in range(slices):
            pdfs_to_merge_lists.append(pdfs_to_merge[start_ind: start_ind+slice_length])
            start_ind += slice_length
        # add the odd rest if something is left:
        if mod != 0:
            pdfs_to_merge_lists.append(pdfs_to_merge[slice_length*slices:])

    # with pdf merger just the pdfs_to_merge list is needed:
    # loop through the slices, prepend already merged pdfs before the single pdfs
    if len(pdfs_to_merge) > slice_length:
        # create all the temporary pdfs slice merge filenames:
        temp_merge_name = [out_path.replace('.pdf', f'_temp_{int(i)}.pdf') for i in range(slices+1)]
        first_pdf = True
        for i in range(len(pdfs_to_merge_lists)):
            merger = PdfMerger()
            first_slice = True
            for pdf in pdfs_to_merge_lists[i]:
                # before any slice, prepend already merged pdfs, unless the first merge
                if first_slice and not first_pdf:
                    first_slice = False
                    first_pdf = False
                    #
                    merger.append(temp_merge_name[i-1])
                    merger.append(pdf)
                else:
                    first_pdf = False
                    first_slice = False
                    merger.append(pdf)
            if i != max(range(len(pdfs_to_merge_lists))):
                merger.write(temp_merge_name[i])
                merger.close()
            else:
                merger.write(out_path)
                merger.close()
        # remove the created temp_merge_name temporary pdf files
        # [os.remove(temp_merge_name[i]) for i in range(len(temp_merge_name)-1) ]
    else:
        merger = PdfMerger()
        for pdf in pdfs_to_merge:
            merger.append(pdf)
        merger.write(out_path)
        merger.close()
