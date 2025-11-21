import type { GraphMakerState } from '@milaboratories/graph-maker';
import type { InferOutputsType, PColumnIdAndSpec, PlDataTableStateV2, PlRef } from '@platforma-sdk/model';
import {
  BlockModel,
  createPlDataTableStateV2,
  createPlDataTableV2,
  isPColumn,
  isPColumnSpec,
} from '@platforma-sdk/model';

export type BlockArgs = {
  datasetRef?: PlRef;
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
    allowRun: false,
  })

  .argsValid((ctx) => {
    const datasetRef = ctx.args.datasetRef;
    return datasetRef !== undefined && datasetRef !== null;
  })

  .output('datasetOptions', (ctx) => {
    return ctx.resultPool.getOptions((v) => {
      if (!isPColumnSpec(v)) return false;

      // Accept all five types: csv/tsv, h5ad, h5, rds (Seurat), and mtx
      const domain = v.domain;
      const hasRoleAxis = v.axesSpec?.some((axis) => axis.name === 'pl7.app/sc/cellRangerFileRole');

      const isCsv = domain !== undefined
        && (domain['pl7.app/fileExtension'] === 'csv'
          || domain['pl7.app/fileExtension'] === 'csv.gz'
          || domain['pl7.app/fileExtension'] === 'tsv'
          || domain['pl7.app/fileExtension'] === 'tsv.gz');

      const isH5ad = domain !== undefined
        && domain['pl7.app/fileExtension'] === 'h5ad';

      const isH5 = domain !== undefined
        && domain['pl7.app/fileExtension'] === 'h5';

      const isRds = domain !== undefined
        && domain['pl7.app/fileExtension'] === 'rds';

      const isMtx = hasRoleAxis;

      return (
        v.name === 'pl7.app/sequencing/data'
        && (v.valueType as string) === 'File'
        && (isCsv || isH5ad || isH5 || isRds || isMtx)
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
