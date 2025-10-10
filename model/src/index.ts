import type { InferOutputsType, PlDataTableStateV2, PlRef } from '@platforma-sdk/model';
import { BlockModel, createPlDataTableStateV2, isPColumnSpec } from '@platforma-sdk/model';

export type BlockArgs = {
  datasetRef?: PlRef;
  name?: string;
};

export type UiState = {
  tableState: PlDataTableStateV2;
  title: string;
  settingsOpen: boolean;
};

export const model = BlockModel.create()

  .withArgs<BlockArgs>({})

  .withUiState<UiState>({
    tableState: createPlDataTableStateV2(),
    title: 'Import GEX Data',
    settingsOpen: true,
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

  .sections((_ctx) => [{ type: 'link', href: '/', label: 'Main' }])

  .title((ctx) => ctx.uiState.title)

  .done(2);

export type BlockOutputs = InferOutputsType<typeof model>;
