import type { GraphMakerState } from '@milaboratories/graph-maker';
import type { InferOutputsType, PColumnIdAndSpec, PlDataTableStateV2, PlRef } from '@platforma-sdk/model';
import { BlockModel, createPlDataTableStateV2, createPlDataTableV2, isPColumn, isPColumnSpec } from '@platforma-sdk/model';

export type BlockArgs = {
  datasetRef?: PlRef;
  mtxDatasetRef?: PlRef;
  h5adDatasetRef?: PlRef;
  importMode: 'csv' | 'mtx' | 'h5ad';
  // name?: string;
};

export type UiState = {
  tableState: PlDataTableStateV2;
  title: string;
  settingsOpen: boolean;
  graphState: GraphMakerState;
  allowRun: boolean;
};

export const model = BlockModel.create()

  .withArgs<BlockArgs>({
    importMode: 'csv',
  })

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
    allowRun: false,
  })

  .argsValid((ctx) => {
    const { importMode, datasetRef, mtxDatasetRef, h5adDatasetRef } = ctx.args;
    if (importMode === 'csv') {
      if (datasetRef === undefined) return false;
      if (!ctx.uiState.allowRun) return false;
    }

    if (importMode === 'mtx') {
      if (mtxDatasetRef === undefined) return false;
    }

    if (importMode === 'h5ad') {
      if (h5adDatasetRef === undefined) return false;
    }

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

  .output('h5adDatasetOptions', (ctx) => {
    return ctx.resultPool.getOptions((v) => {
      if (!isPColumnSpec(v)) return false;
      const domain = v.domain;
      return (
        v.name === 'pl7.app/sequencing/data'
        && (v.valueType as string) === 'File'
        && domain !== undefined
        && (domain['pl7.app/fileExtension'] === 'h5ad')
      );
    },
    );
  })

  .output('mtxDatasetOptions', (ctx) => {
    return ctx.resultPool.getOptions((v) => {
      if (!isPColumnSpec(v)) return false;
      // Check if it has the cellRangerFileRole axis
      const hasRoleAxis = v.axesSpec?.some((axis) => axis.name === 'pl7.app/sc/cellRangerFileRole');
      return (
        v.name === 'pl7.app/sequencing/data'
        && (v.valueType as string) === 'File'
        && hasRoleAxis
      );
    },
    );
  })

  .output('errorLog', (ctx) => {
    const pCols = ctx.prerun?.resolve({ field: 'errorLog', allowPermanentAbsence: true })?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    return ctx.createPFrame(pCols);
  })

  .output('resultsSummaryPf', (ctx) => {
    const pCols = ctx.outputs?.resolve('resultsSummaryPf')?.getPColumns();
    if (pCols === undefined) {
      return undefined;
    }

    return createPlDataTableV2(ctx, pCols, ctx.uiState.tableState);
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

  .output('isRunning', (ctx) => ctx.outputs?.getIsReadyOrError() === false
    || ctx.prerun?.getIsReadyOrError() === false)

  .output('runningPrerun', (ctx) => ctx.prerun?.getIsReadyOrError() === false)

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
