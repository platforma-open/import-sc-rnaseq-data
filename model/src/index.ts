import type { InferOutputsType, PColumnIdAndSpec, PlDataTableStateV2, PlRef } from '@platforma-sdk/model';
import { BlockModel, createPlDataTableStateV2, isPColumn, isPColumnSpec } from '@platforma-sdk/model';
import type { GraphMakerState } from '@milaboratories/graph-maker';

export type BlockArgs = {
  datasetRef?: PlRef;
  name?: string;
};

export type UiState = {
  tableState: PlDataTableStateV2;
  title: string;
  settingsOpen: boolean;
  graphState: GraphMakerState;
};

export const model = BlockModel.create()

  .withArgs<BlockArgs>({})

  .withUiState<UiState>({
    tableState: createPlDataTableStateV2(),
    title: 'Import GEX Data',
    settingsOpen: true,
    graphState: {
      template: 'violin',
      title: 'Cell QC metrics',
      // layersSettings: {
      //   violin: {
      //     fillColor: '#99E099',
      //   },
      // },
    },
  })

  .argsValid((ctx) => {
    const { datasetRef } = ctx.args;
    if (datasetRef === undefined) return false;

    return true;
  })

  .output('datasetOptions', (ctx) => {
    return ctx.resultPool.getOptions((v) => {
      if (!isPColumnSpec(v)) return false;
      const domain = v.domain;
      return (
        v.name === 'pl7.app/sequencing/data'
        && (v.valueType as string) === 'File'
        && domain !== undefined
        && (domain['pl7.app/fileExtension'] === 'csv'
          || domain['pl7.app/fileExtension'] === 'csv.gz'
          || domain['pl7.app/fileExtension'] === 'tsv'
          || domain['pl7.app/fileExtension'] === 'tsv.gz')
      );
    },
    );
  })

  .output('cellMetricsPf', (wf) => {
    const pCols = wf.outputs?.resolve('cellMetricsPf')?.getPColumns();
    if (pCols === undefined) return undefined;

    const upstream = wf.resultPool
      .getData()
      .entries.map((v) => v.obj)
      .filter(isPColumn)
      .filter((col) => col.spec.name === 'pl7.app/label');

    return wf.createPFrame([...pCols, ...upstream]);
  })

  .output('cellMetricsPfDefaults', (wf) => {
    let pCols = wf.outputs?.resolve('cellMetricsPf')?.getPColumns();
    if (pCols === undefined) return undefined;

    // Add sample labels
    const upstream = wf.resultPool
      .getData()
      .entries.map((v) => v.obj)
      .filter(isPColumn)
      .filter((col) => col.spec.name === 'pl7.app/label');

    pCols = [...pCols, ...upstream];
    return pCols.map(
      (c) =>
        ({
          columnId: c.id,
          spec: c.spec,
        } satisfies PColumnIdAndSpec),
    );
  })

  .output('cellMetricsSpec', (wf) => {
    const pCols = wf.outputs?.resolve('cellMetricsPf')?.getPColumns();
    if (pCols === undefined) return undefined;
    return pCols[0].spec;
  })

  .output('isRunning', (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .sections([
    { type: 'link', href: '/', label: 'Main' },
    { type: 'link', href: '/CellQC', label: 'Cell QC' },
  ])

  .title((ctx) =>
    ctx.uiState.title
      ? `Import scRNA-seq Data - ${ctx.uiState.title}`
      : 'Import scRNA-seq Data',
  )

  .done(2);

export type BlockOutputs = InferOutputsType<typeof model>;
