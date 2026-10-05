# LDC post-processing

Pasta isolada para gerar figuras combinadas de runs LDC sem mexer nos plots em
`postprocess/`.

Por padrao os scripts leem:

- `postprocess_ldc/ids.txt`
- `runs/<run_id>/outputs/centerline.bin`
- `runs/<run_id>/outputs/tke.bin`
- `runs/<run_id>/vtk/output_*.vti`

O formato do `ids.txt` aceita um run por linha. Tambem aceita label opcional:

```text
20260720_174720_D2Q9_double_256x64_RE100_DY0p0_T200_GPU0 | Couette 256x64
```

## Rodar tudo

```bash
python3 postprocess_ldc/plot_all.py
```

## Rodar partes separadas

```bash
python3 postprocess_ldc/plot_centerlines.py
python3 postprocess_ldc/plot_tke.py
python3 postprocess_ldc/plot_vti_velocity.py
```

As figuras saem em `postprocess_ldc/figures/`.

## Opcoes uteis

Usar outro arquivo de ids:

```bash
python3 postprocess_ldc/plot_all.py --ids caminho/ids.txt
```

Plotar todos os VTIs, em vez de apenas o ultimo de cada run:

```bash
python3 postprocess_ldc/plot_vti_velocity.py --which all
```

Gerar campo com streamlines alem de magnitude:

```bash
python3 postprocess_ldc/plot_vti_velocity.py --streamlines
```

## Dependencias

```bash
python3 -m pip install -r postprocess_ldc/requirements.txt
```
