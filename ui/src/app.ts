import { model } from '@platforma-open/milaboratories.import-sc-rnaseq-data.model';
import { defineApp } from '@platforma-sdk/ui-vue';
import MainPage from './pages/MainPage.vue';
import CellQC from './pages/CellQC.vue';

export const sdkPlugin = defineApp(model, (app) => {
  return {
    progress: () => {
      return app.model.outputs.isRunning;
    },
    showErrorsNotification: true,
    routes: {
      '/': () => MainPage,
      '/CellQC': () => CellQC,
    },
  };
});

export const useApp = sdkPlugin.useApp;
